using Test

using DiffEqBase, ForwardDiff
using DiffEqBase: prob2dtmin, timedepentdtmin, _rate_prototype
using Unitful
using ForwardDiff: Dual, Tag

@testset "tspan2dtmin" begin
    # we only need to test very rough equality since timestepping isn't science.
    function approxoftype(a, b; rtol = 0.5)
        return typeof(a) === typeof(b) && isapprox(a, b; rtol = rtol)
    end
    function tspan2dtmin(tspan; kwargs...)
        prob2dtmin(ODEProblem((u, p, t) -> u, 1, tspan); kwargs...)
    end
    @test approxoftype(tspan2dtmin((10, 100.0)), eps(100.0))
    @test approxoftype(tspan2dtmin((-10000.0, 100.0)), eps(10000.0))
    @test tspan2dtmin((1, 2)) === 0
    @test approxoftype(tspan2dtmin((1 // 10, 2 // 10)), 1 // 2^33)
    @test approxoftype(tspan2dtmin((2 // 10, Inf)), eps(1.0))
    @test approxoftype(tspan2dtmin((2 // 1, Inf)), eps(2.0))
    @test approxoftype(tspan2dtmin((0, Inf)), eps(1.0))
    @test approxoftype(tspan2dtmin((0.0, Inf)), eps(1.0))
    @test approxoftype(tspan2dtmin((0.0, 1.0e-6)), eps(1.0e-6))
    @test approxoftype(tspan2dtmin((1.0e6, 1.0e6 + 1)), eps(1.0e6))
    @test_throws ArgumentError tspan2dtmin((Inf, 100.0))
    @test approxoftype(tspan2dtmin((0.0f0, 1.0f5); use_end_time = false), eps(1.0f0))
    @test approxoftype(timedepentdtmin(10.0f0, eps(1.0f0)), eps(10.0f0))
    @test approxoftype(timedepentdtmin(10, eps(1.0f0)), eps(1.0f0))
end

@testset "prob2dtmin" begin
    @test prob2dtmin((0.0, 10.0), 1.0, false) == eps(Float64)
    @test prob2dtmin((0.0f0, 10.0f0), 1.0f0, false) == eps(Float32)
    @test prob2dtmin((0.0, 10.0), ForwardDiff.Dual(1.0), false) == eps(Float64)
end

@testset "_rate_prototype" begin
    @test _rate_prototype([1.0f0], 1.0, 1.0) isa Vector{Float32}
    td = Dual{Tag{typeof(+), Float64}}(2.0, 1.0)
    @test _rate_prototype([1.0f0], td, td) isa Vector{Float32}
    xd = [Dual{Tag{typeof(+), Float32}}(2.0, 1.0)]
    @test _rate_prototype(xd, 1.0, 1.0) isa typeof(xd)
    @test _rate_prototype([u"1f0m"], u"1.0s", 1.0) isa typeof([u"1f0m/s"])
end
