# Regression test for the DiffEqBaseEnzymeExt `solve_up` reverse rule guard
# (SciML/OrdinaryDiffEq.jl#3740): under `set_runtime_activity`, a u0 that aliases a
# Const-annotated problem's own array must not receive the du0 cotangent (its
# "shadow" IS the primal). First gradient call used to be correct while silently
# corrupting `prob.u0`; the second call then returned garbage.
#
# Deps: SciMLSensitivity (real adjoint), OrdinaryDiffEqVerner, Enzyme, ForwardDiff,
# ModelingToolkit, ModelingToolkitBase, SymbolicIndexingInterface.

using DiffEqBase, Enzyme, Test
using SciMLSensitivity
using OrdinaryDiffEqVerner
using ForwardDiff

f_oop(u, p, t) = p .* u
u0_init = [2.0, 3.0]
p_init = [0.5, 0.7]
prob = ODEProblem(f_oop, copy(u0_init), (0.0, 1.0), copy(p_init))

solve_kwargs = (; saveat = 0.25, abstol = 1.0e-8, reltol = 1.0e-8)

@testset "runtime-activity aliased u0 is not corrupted" begin
    # u0 is the Const problem's own array reused via remake — the
    # `setsym_oop`/`remake` pattern of MTK loss functions.
    function loss_aliased(p, q)
        prob = q[1]
        prob2 = remake(prob; prob.u0, p)
        sol = solve(prob2, Vern7(); sensealg = GaussAdjoint(), solve_kwargs...)
        return sum(abs2, Array(sol))
    end

    q = (prob,)
    g_ref = ForwardDiff.gradient(p -> loss_aliased(p, q), p_init)

    g1 = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), loss_aliased, copy(p_init), Const(q))[1]
    @test g1 ≈ g_ref rtol = 1.0e-5
    @test prob.u0 == u0_init     # primal problem must NOT have been mutated
    g2 = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), loss_aliased, copy(p_init), Const(q))[1]
    @test g2 ≈ g_ref rtol = 1.0e-5 # second call sees uncorrupted state
end

@testset "genuinely active u0 still accumulates" begin
    # The guard must not skip accumulation when a real shadow exists: here both
    # u0 and p derive from the differentiated input, so du0 must flow.
    function loss_active(x, q)
        prob2 = remake(q[1]; u0 = x[1:2], p = x[3:4])
        sol = solve(prob2, Vern7(); sensealg = GaussAdjoint(), solve_kwargs...)
        return sum(abs2, Array(sol))
    end

    q = (prob,)
    x0 = vcat(u0_init, p_init)
    g_ref = ForwardDiff.gradient(x -> loss_active(x, q), x0)
    @test maximum(abs, g_ref[1:2]) > 0  # sanity: u0 gradient is nonzero

    g = Enzyme.gradient(set_runtime_activity(Enzyme.Reverse), loss_active, copy(x0), Const(q))[1]
    @test g ≈ g_ref rtol = 1.0e-5
end

@testset "parameter derivative read from the solution's problem" begin
    # A loss can depend on `p` through the returned solution, e.g. an observed variable
    # evaluated with `sol.prob.p` (SciML/SciMLSensitivity.jl#1714). Enzyme accumulates
    # that part into the shadow of `sol.prob.p`, which the rule must route back to `p`
    # exactly once.
    function loss_prob_p(p, q, sensealg)
        sol = solve(remake(q[1]; p), Vern7(); sensealg, solve_kwargs...)
        return sum(abs2, Array(sol)) + sum(abs2, sol.prob.p)
    end

    q = (prob,)
    @testset "$(nameof(typeof(sensealg)))" for sensealg in (
            GaussAdjoint(autojacvec = EnzymeVJP()),
            InterpolatingAdjoint(autojacvec = EnzymeVJP()),
            QuadratureAdjoint(autojacvec = EnzymeVJP()),
        )
        g_ref = ForwardDiff.gradient(p -> loss_prob_p(p, q, nothing), p_init)
        g1 = Enzyme.gradient(
            set_runtime_activity(Enzyme.Reverse), loss_prob_p, copy(p_init), Const(q),
            Const(sensealg)
        )[1]
        @test g1 ≈ g_ref rtol = 1.0e-5
        g2 = Enzyme.gradient(
            set_runtime_activity(Enzyme.Reverse), loss_prob_p, copy(p_init), Const(q),
            Const(sensealg)
        )[1]
        @test g2 ≈ g1
    end
end

# The case from SciML/SciMLSensitivity.jl#1714: an MTK observed variable read from the
# solution. MTK problems default to `AutoDespecialize`, so the shadow of `sol.prob.p` is a
# `DespecializedParameters` around an `MTKParameters`, and the rule has to route its
# tunable buffer back to `p`.
using ModelingToolkit, SymbolicIndexingInterface
import SciMLBase
using ModelingToolkitBase: MTKParameters
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters k = 1.0 c = 0.5
@variables x(t) = 1.0 y(t)
obs_sys = mtkcompile(System([D(x) ~ -k * x, y ~ x + c], t; name = :obs))
obs_prob = ODEProblem(obs_sys, [], (0.0, 2.0))
const obs_setter = setsym_oop(obs_prob, [obs_sys.k, obs_sys.c])
const get_x = getsym(obs_prob, obs_sys.x)
const get_y = getsym(obs_prob, obs_sys.y)

function obs_solve(θ, prob, sensealg)
    u0, p = obs_setter(prob, θ)
    return solve(remake(prob; u0, p), Vern7(); saveat = 0.1, sensealg, solve_kwargs...)
end
# `y = x + c` over 21 save points, so ∂L/∂c = 21 exactly; `c` only enters through `y`.
loss_observed(θ, prob, sensealg) = sum(get_y(obs_solve(θ, prob, sensealg)))
# Mixed: the state part goes through the adjoint, the observed part through `sol.prob.p`.
function loss_mixed(θ, prob, sensealg)
    sol = obs_solve(θ, prob, sensealg)
    return sum(abs2, get_x(sol)) + sum(abs2, get_y(sol))
end

@testset "MTK observed variable read from the solution" begin
    @test SciMLBase.unwrap_parameters(obs_prob.p) isa MTKParameters
    θ = [1.0, 0.5]
    @testset "$(sensealg === nothing ? "automatic" : nameof(typeof(sensealg)))" for sensealg in (
            nothing,
            GaussAdjoint(autojacvec = EnzymeVJP()),
            InterpolatingAdjoint(autojacvec = EnzymeVJP()),
            QuadratureAdjoint(autojacvec = EnzymeVJP()),
        )
        @testset "$(nameof(loss))" for loss in (loss_observed, loss_mixed)
            g_ref = ForwardDiff.gradient(θ -> loss(θ, obs_prob, nothing), θ)
            g1 = Enzyme.gradient(
                set_runtime_activity(Enzyme.Reverse), loss, copy(θ), Const(obs_prob),
                Const(sensealg)
            )[1]
            @test g1 ≈ g_ref rtol = 1.0e-5
            loss === loss_observed && @test g1[2] ≈ 21 rtol = 1.0e-8
            g2 = Enzyme.gradient(
                set_runtime_activity(Enzyme.Reverse), loss, copy(θ), Const(obs_prob),
                Const(sensealg)
            )[1]
            @test g2 ≈ g1
        end
    end
end
