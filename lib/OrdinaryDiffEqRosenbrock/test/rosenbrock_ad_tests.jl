using OrdinaryDiffEqRosenbrock, Test
using SciMLBase: ODEProblem
using ForwardDiff

# Regression test for issue #3486 — nested ForwardDiff (hessian /
# gradient-of-gradient) through a Rosenbrock solve used to error with:
#
#   MethodError: no method matching Float64(::ForwardDiff.Dual{...})
#     @ setproperty!(::JacReuseState{Float64}, ::Symbol, ::Dual)
#
# The Jacobian-reuse state was allocated via
# `zero(constvalue(uBottomEltypeNoUnits))`, which unconditionally strips
# ForwardDiff.Dual layers and pins the dtgamma field type to Float64. Under
# `ForwardDiff.hessian`, dtgamma is a nested `Dual{Outer,Dual{Inner,Float64}}`
# and can't be stored into a Float64 field. The fix seeds the state from
# `zero(dt)` instead, so the field type matches dtgamma's true solve-time type.
#
# Cover both OOP and IIP code paths, and both styles of constant/mutable
# cache:
#   * Rosenbrock23 / Rosenbrock32: their own constant- and mutable-cache types
#     (Rosenbrock23ConstantCache / Rosenbrock23Cache etc.).
#   * RodasTableauAlgorithms: shared RosenbrockCombinedConstantCache (OOP) and
#     RosenbrockCache (IIP).
# Spot-check across several RodasTableau methods so different stage counts /
# interp_order branches get exercised.

# Analytic reference for u' = -p₁*u², u(0)=1, integrated to t=p₂:
#   u(t) = 1 / (1 + p₁*t)
#   ∂u/∂p₁ = -t / (1 + p₁*t)²                               # at p₁=p₂=1 → -1/4
#   ∂u/∂p₂ = -p₁ / (1 + p₁*t)²                              # at p₁=p₂=1 → -1/4
#   ∂²u/∂p₁² = 2t² / (1 + p₁*t)³                            # at p₁=p₂=1 →  1/4
#   ∂²u/∂p₂² = 2p₁² / (1 + p₁*t)³                           # at p₁=p₂=1 →  1/4
#   ∂²u/∂p₁∂p₂ = (2p₁*t - 1 - p₁*t) / (1 + p₁*t)³ (mixed)   # at p₁=p₂=1 →  0
# So at p = [1,1], gradient ≈ [-1/4, -1/4] and hessian ≈ [1/4 0; 0 1/4].
rhs_oop_simple(u, p, t) = [-p[1] * u[1] * u[1]]
function rhs_iip_simple!(du, u, p, t)
    du[1] = -p[1] * u[1] * u[1]
    return nothing
end

# Tight tolerances so the lower-order methods (Rosenbrock23 etc.) land close
# enough to the analytic values that this test reads as "did AD actually
# propagate through correctly" rather than "did the integrator's default
# tolerance happen to be tight enough".
function make_loss(alg, ::Val{iip}) where {iip}
    if iip
        # u0 must be promoted to eltype(p) so the in-place writes don't blow up
        # on a Float64 buffer when p (and hence du) is Dual.
        return p -> begin
            u0 = [one(eltype(p))]
            prob = ODEProblem{true}(rhs_iip_simple!, u0, (zero(eltype(p)), p[2]), p)
            sol = solve(prob, alg; abstol = 1.0e-10, reltol = 1.0e-10)
            only(last(sol.u))
        end
    else
        return p -> begin
            prob = ODEProblem{false}(rhs_oop_simple, [1.0], (0.0, p[2]), p)
            sol = solve(prob, alg; abstol = 1.0e-10, reltol = 1.0e-10)
            only(last(sol.u))
        end
    end
end

const ROSENBROCK_ALGS_TO_TEST = (
    # Rosenbrock23 / Rosenbrock32 each have their own ConstantCache type.
    Rosenbrock23(),
    Rosenbrock32(),
    # RodasTableauAlgorithms all share RosenbrockCombinedConstantCache —
    # sample a spread (2nd-, 3rd-, 4th-, 5th-order; W-method vs strict;
    # stiffly-accurate vs not). Leave out the ROSnnPWX family here because
    # several of them hit dtmin on this simple non-stiff problem under
    # nested ForwardDiff and that failure mode is not what we're testing.
    Rodas3(),
    Rodas23W(),
    ROS3P(),
    Rodas4(),
    Rodas5P(),
)

@testset "Nested ForwardDiff through solve (issue #3486) — $(nameof(typeof(alg))) iip=$iip" for alg in
                                                                                                                                                        ROSENBROCK_ALGS_TO_TEST,
        iip in (false, true)

    loss = make_loss(alg, Val(iip))

    # First-order differentiation already worked before the fix — included as
    # a correctness anchor.
    g = ForwardDiff.gradient(loss, [1.0, 1.0])
    @test all(isfinite, g)
    @test g ≈ [-0.25, -0.25] atol = 1.0e-4

    # Second-order differentiation is the one that used to throw.
    H = ForwardDiff.hessian(loss, [1.0, 1.0])
    @test all(isfinite, H)
    @test H ≈ [0.25 0.0; 0.0 0.25] atol = 1.0e-4
end
