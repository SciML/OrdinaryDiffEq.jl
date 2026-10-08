using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, JLArrays, StaticArrays, Test
using RecursiveArrayTools: ArrayPartition
using ComponentArrays: ComponentArray
using ArrayInterface: ArrayInterface

f(u, p, t) = @. -u / 2
build_vcc(affect!) = VectorContinuousCallback(
    (out, u, t, integrator) -> (out .= u .- 0.5),
    affect!, 2
)

@testset "Out-of-place VectorContinuousCallback with JLArray u" begin
    JLArrays.allowscalar(false)
    u0 = JLArray([1.0, 2.0])
    sol = solve(ODEProblem(f, u0, (0.0, 1.0)), Tsit5(); callback = build_vcc((integrator, idx) -> nothing))
    @test sol.retcode == ReturnCode.Success
    @test sol.u[end] isa JLArray
end

@testset "Vector u gives the analytic event time" begin
    # u0 = [1.0, 2.0] decays as exp(-t/2); component 1 crosses 0.5 at t = 2log(2).
    u0 = [1.0, 2.0]
    sol = solve(ODEProblem(f, u0, (0.0, 5.0)), Tsit5(); callback = build_vcc((integrator, idx) -> terminate!(integrator)))
    @test sol.retcode == ReturnCode.Terminated
    @test sol.t[end] ≈ 2 * log(2) rtol = 1.0e-4
end

@testset "_build_callback_cache tracks ArrayInterface.ismutable, not Base.ismutable" begin
    # ArrayPartition and a strided view of a JLArray are themselves immutable wrapper
    # structs (`Base.ismutable == false`) even though their underlying storage can be
    # mutated (`ArrayInterface.ismutable == true`). The out-of-place cache must use
    # `ArrayInterface.ismutable` so a JLArray wrapped in either still gets a cache
    # built from `similar(u, ...)` (a JLArray, or an ArrayPartition of JLArrays)
    # instead of a plain CPU `Array` from `zeros`.
    for u in (
            ArrayPartition(JLArray([1.0]), JLArray([2.0])),
            view(JLArray(collect(1.0:4.0)), 1:2:3),
        )
        @test !Base.ismutable(u)
        @test ArrayInterface.ismutable(u)
        cache = OrdinaryDiffEqCore._build_callback_cache(u, 2, Val(false), Float64)
        @test !(cache.tmp_condition isa Array)
    end
end

@testset "Out-of-place VectorContinuousCallback with a strided SubArray of JLArray" begin
    JLArrays.allowscalar(false)
    u0 = view(JLArray(collect(1.0:4.0)), 1:2:3)
    sol = solve(ODEProblem(f, u0, (0.0, 1.0)), Tsit5(); callback = build_vcc((integrator, idx) -> nothing))
    @test sol.retcode == ReturnCode.Success
end

@testset "Out-of-place VectorContinuousCallback with ComponentArray u" begin
    u0 = ComponentArray(a = 1.0, b = 2.0)
    sol = solve(ODEProblem(f, u0, (0.0, 1.0)), Tsit5(); callback = build_vcc((integrator, idx) -> nothing))
    @test sol.retcode == ReturnCode.Success
end

@testset "Out-of-place VectorContinuousCallback with SVector u is unaffected" begin
    # SVector is immutable under both `Base.ismutable` and `ArrayInterface.ismutable`,
    # so `_build_callback_cache` keeps using the CPU-`zeros` constructor unchanged.
    u0 = SVector(1.0, 2.0)
    sol = solve(ODEProblem(f, u0, (0.0, 1.0)), Tsit5(); callback = build_vcc((integrator, idx) -> nothing))
    @test sol.retcode == ReturnCode.Success
    @test sol.u[end] ≈ u0 .* exp(-1.0 / 2) rtol = 1.0e-6
end
