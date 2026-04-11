using Test
using OrdinaryDiffEqRosenbrock
using SciMLBase
using ADTypes
using ForwardDiff

# Regression test for nested ForwardDiff over an ODE solve (#3381). When a
# Rosenbrock solver is invoked with a `Vector{<:Dual}` `p` (i.e. we are
# inside some *outer* ForwardDiff layer), the inner Rosenbrock Jacobian
# widens `u` to a deeper nested-Dual type via its JacobianConfig. The user
# `f(du, u, p, t)` body then multiplies `p[i] * u[i]` across two different
# Dual nesting levels. ForwardDiff's cross-tag multiplication uses
# `tagcount`-based precedence, which is a `@generated` function whose
# literal value is baked at first compile. Depending on the order in which
# tag types are first instantiated (which varies with precompile order),
# the resulting ordering can invert nesting and produce a triple-nested
# `Dual` that crashes the eventual `setindex!` with
#     `Float64(::nested_dual)` MethodError.
#
# The fix (`_widen_uf_p_for_jac` in derivative_wrappers.jl) lifts
# `uf.p` into the inner nested-Dual type ahead of time so the user
# body never has to cross tag levels.
#
# To reliably reproduce the broken precedence without depending on
# NonlinearSolve's precompile-baked NonlinearSolveTag, we construct the
# outer Dual manually with a throwaway tag that is *not* a concrete
# `Tag(f, V)` call. Because the throwaway tag is only used as a type
# parameter, ForwardDiff does not trigger its `tagcount` until the first
# cross-tag multiplication — which happens *after* `OrdinaryDiffEqTag` has
# already had its nested-V `tagcount` triggered via
# `OrdinaryDiffEqDifferentiation.prepare_ADType`. This matches the
# precedence inversion that the original issue hits via
# `NonlinearSolveTag`'s precompile workload.

struct _Nested3381OuterTag end

function _nested3381_ode!(du, u, p, t)
    du[1] = -p[1] * u[1]
    du[2] = -u[1] - p[2] * u[2]
    return nothing
end

@testset "Nested ForwardDiff through Rosenbrock Jacobian (#3381)" begin
    OuterTag = ForwardDiff.Tag{_Nested3381OuterTag, Float64}
    DualT = ForwardDiff.Dual{OuterTag, Float64, 2}
    pdual = [
        DualT(1.5, ForwardDiff.Partials((1.0, 0.0))),
        DualT(2.0, ForwardDiff.Partials((0.0, 1.0))),
    ]

    ode_f = ODEFunction{true, SciMLBase.FullSpecialize}(_nested3381_ode!)
    u0 = [1.0, 1.0]
    tspan = (0.0, 1.0)
    prob = ODEProblem(ode_f, u0, tspan, pdual)

    sol = solve(prob, Rosenbrock23(autodiff = AutoForwardDiff(chunksize = 2)))
    @test SciMLBase.successful_retcode(sol.retcode)
    @test eltype(sol.u[end]) === DualT
    # Sanity: primal values still reach the target trajectory.
    @test all(isfinite, ForwardDiff.value.(sol.u[end]))
end
