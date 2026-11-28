using OrdinaryDiffEqDifferentiation: dae_jvp_cache, set_dae_jvp_cj!
using SciMLOperators: update_coefficients!
using SciMLBase: DAEFunction
using ADTypes
using ADTypes: AutoForwardDiff, AutoFiniteDiff
using LinearAlgebra: mul!
using Test

# Nonlinear in `du`. At `u = [1, 1]`, `du = [-1, -2]` the Jacobian of the residual in `u`
# plus `cj` times its Jacobian in `du` is `[3(cj + 1) 0; -2 1]`.
function manufactured!(r, d, u, p, t)
    r[1] = (d[1] + u[1]) * (1 + d[1]^2 + u[2]^2)
    r[2] = u[2] - u[1]^2
    return nothing
end

@testset "Matrix-free DAE Newton operator at a large cj" begin
    u, du, v = [1.0, 1.0], [-1.0, -2.0], [1.0, 0.5]
    for autodiff in (AutoForwardDiff(), AutoFiniteDiff(; fdtype = Val(:central))), cj in (1.0, 1.0e6, 1.0e10)
        J = dae_jvp_cache(DAEFunction{true}(manufactured!), copy(u), nothing, 0.0; autodiff)
        update_coefficients!(J, copy(u), nothing, 0.0; du = copy(du))
        set_dae_jvp_cj!(J, cj)
        Jv = mul!(zeros(2), J, v)
        # The algebraic component is `-1.5` however large `cj` is, so it is held to its own
        # tolerance rather than to one scaled by the differential component.
        @test Jv[1] ≈ 3 * (cj + 1) rtol = 1.0e-6
        @test Jv[2] ≈ -1.5 rtol = 1.0e-6
    end
end

@testset "Matrix-free DAE Newton operator rejects non-centered finite differences" begin
    f = DAEFunction{true}(manufactured!)
    for autodiff in (
            AutoFiniteDiff(), AutoFiniteDiff(; fdtype = Val(:complex)),
            ADTypes.AutoSparse(AutoFiniteDiff()), ADTypes.AutoFiniteDifferences(; fdm = nothing),
        )
        @test_throws ArgumentError dae_jvp_cache(f, [1.0, 1.0], nothing, 0.0; autodiff)
    end
    centered = ADTypes.AutoSparse(AutoFiniteDiff(; fdtype = Val(:central)))
    @test dae_jvp_cache(f, [1.0, 1.0], nothing, 0.0; autodiff = centered) !== nothing
end
