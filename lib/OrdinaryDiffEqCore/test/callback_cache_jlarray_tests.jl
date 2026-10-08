using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, JLArrays, StaticArrays, Test

@testset "Out-of-place VectorContinuousCallback with JLArray u" begin
    JLArrays.allowscalar(false)
    f(u, p, t) = @. -u / 2
    u0 = JLArray([1.0, 2.0])
    vcc = VectorContinuousCallback(
        (out, u, t, integrator) -> (out .= u .- 0.5),
        (integrator, idx) -> nothing, 2
    )
    sol = solve(ODEProblem(f, u0, (0.0, 1.0)), Tsit5(); callback = vcc)
    @test sol.retcode == ReturnCode.Success
    @test sol.u[end] isa JLArray
end

@testset "Out-of-place VectorContinuousCallback with SVector u is unaffected" begin
    # SVector u is immutable, so `_build_callback_cache` keeps using the
    # pre-existing `DiffEqBase.CallbackCache(max_len, T, T)` (CPU `zeros`) path
    # regardless of this fix; this just guards that that path still works.
    f(u, p, t) = @. -u / 2
    u0 = SVector(1.0, 2.0)
    vcc = VectorContinuousCallback(
        (out, u, t, integrator) -> (out .= u .- 0.5),
        (integrator, idx) -> nothing, 2
    )
    sol = solve(ODEProblem(f, u0, (0.0, 1.0)), Tsit5(); callback = vcc)
    @test sol.retcode == ReturnCode.Success
    @test sol.u[end] ≈ u0 .* exp(-1.0 / 2) rtol = 1.0e-6
end
