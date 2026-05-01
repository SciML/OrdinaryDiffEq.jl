using Test
using OrdinaryDiffEqRosenbrock
using SciMLBase
using ADTypes
using ForwardDiff

# Regression test for nested ForwardDiff over an ODE solve
# (https://github.com/SciML/OrdinaryDiffEq.jl/issues/3381).
#
# When a Rosenbrock solver is invoked with a `Vector{<:Dual}` `p` (i.e. we are
# inside an *outer* ForwardDiff layer), the inner Rosenbrock Jacobian widens
# `u` into a deeper nested-Dual type via its `JacobianConfig`. The user body
# `f(du, u, p, t)` then multiplies `p[i] * u[i]` across two different Dual
# nesting levels. ForwardDiff's cross-tag multiplication goes through a
# `@generated tagcount` precedence whose literal value is baked at first
# compile and depends on precompile ordering; that ordering can invert the
# nesting and produce a triple-nested `Dual` that crashes the eventual
# `setindex!(du, ...)` with `Float64(::nested_dual)`.
#
# The fix (`_widen_uf_p_for_jac` in derivative_wrappers.jl) lifts `uf.p` into
# the inner nested-Dual type ahead of `DI.jacobian!`, so the user body never
# multiplies across tag levels.
@testset "Nested ForwardDiff through Rosenbrock (issue #3381)" begin
    function ode!(du, u, p, t)
        du[1] = -p[1] * u[1]
        du[2] = -u[1] - p[2] * u[2]
        return nothing
    end
    ode_f = ODEFunction{true, SciMLBase.FullSpecialize}(ode!)

    outer_f = function (p)
        prob = ODEProblem(ode_f, [1.0, 1.0], (0.0, 1.0), p)
        sol = solve(
            prob, Rosenbrock23(autodiff = AutoForwardDiff(chunksize = 2));
            reltol = 1.0e-8, abstol = 1.0e-8
        )
        return sol.u[end]
    end

    J = ForwardDiff.jacobian(outer_f, [1.5, 2.0])
    @test size(J) == (2, 2)
    @test all(isfinite, J)

    # Nested case with a hand-rolled outer tag — mimics NonlinearSolve's
    # Tag{NonlinearSolveTag, Float64} wrapping `p` while Rosenbrock seeds `u`
    # under Tag{OrdinaryDiffEqTag, …}.
    T = ForwardDiff.Tag{:NestedForwardDiffOuter, Float64}
    p_dual = ForwardDiff.Dual{T, Float64, 2}[
        ForwardDiff.Dual{T}(1.5, ForwardDiff.Partials{2, Float64}((1.0, 0.0))),
        ForwardDiff.Dual{T}(2.0, ForwardDiff.Partials{2, Float64}((0.0, 1.0))),
    ]
    prob = ODEProblem(ode_f, [1.0, 1.0], (0.0, 1.0), p_dual)
    sol = solve(
        prob, Rosenbrock23(autodiff = AutoForwardDiff(chunksize = 2));
        reltol = 1.0e-8, abstol = 1.0e-8
    )
    @test SciMLBase.successful_retcode(sol)
    @test all(u -> all(isfinite, u), sol.u)

    # AutoSpecialize wraps `ode!` in a FunctionWrappersWrapper. DiffEqBase's
    # `wrapfun_iip` compiles a Jacobian-case slot with the promoted `p` type;
    # the widen step in `jacobian!` then dispatches into that slot whose body
    # multiplies `u*p` within one tag hierarchy.
    ode_f_auto = ODEFunction{true, SciMLBase.AutoSpecialize}(ode!)
    prob_auto = ODEProblem(ode_f_auto, [1.0, 1.0], (0.0, 1.0), p_dual)
    sol_auto = solve(
        prob_auto, Rosenbrock23(autodiff = AutoForwardDiff(chunksize = 2));
        reltol = 1.0e-8, abstol = 1.0e-8
    )
    @test SciMLBase.successful_retcode(sol_auto)
    @test all(u -> all(isfinite, u), sol_auto.u)
end
