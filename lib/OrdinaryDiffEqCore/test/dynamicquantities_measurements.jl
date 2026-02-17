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

    # Known-broken regression (should become a normal @test when the underlying issue is fixed).
    ok = true
    try
        sol = solve(prob, Tsit5(); abstol = 1e-9, reltol = 1e-9)

        ok &= sol.u[end] isa typeof(u0)
        ok &= eltype(sol.u) == typeof(u0)

        uend_m = ustrip(u"m", sol.u[end])
        ok &= isapprox(Measurements.value(uend_m), exp(1); atol = 1e-6)
        ok &= Measurements.uncertainty(uend_m) > 0
    catch
        ok = false
    end

    @test_broken ok
end
