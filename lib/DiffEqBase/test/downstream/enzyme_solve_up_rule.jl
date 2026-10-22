# Regression test for the DiffEqBaseEnzymeExt `solve_up` reverse rule guard
# (SciML/OrdinaryDiffEq.jl#3740): under `set_runtime_activity`, a u0 that aliases a
# Const-annotated problem's own array must not receive the du0 cotangent (its
# "shadow" IS the primal). First gradient call used to be correct while silently
# corrupting `prob.u0`; the second call then returned garbage.
#
# Deps: SciMLSensitivity (real adjoint), OrdinaryDiffEqVerner, Enzyme, ForwardDiff.

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
