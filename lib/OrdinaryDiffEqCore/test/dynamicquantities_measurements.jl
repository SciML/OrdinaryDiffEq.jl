using DiffEqBase
using OrdinaryDiffEqCore
using OrdinaryDiffEqTsit5: Tsit5
using DynamicQuantities
using Measurements
using Test

@testset "DynamicQuantities units + Measurements uncertainty" begin
    u0 = (1.0 ± 0.1) * (1.0u"m")
    tspan = (0.0u"s", 1.0u"s")

    f(u, p, t) = u / (1u"s")
    prob = ODEProblem(f, u0, tspan)

    sol = solve(prob, Tsit5(); abstol = 1e-9, reltol = 1e-9)

    @test sol.u[end] isa typeof(u0)
    @test eltype(sol.u) == typeof(u0)

    uend_m = ustrip(u"m", sol.u[end])
    @test Measurements.value(uend_m) ≈ exp(1) atol = 1e-6
    @test Measurements.uncertainty(uend_m) > 0
end
