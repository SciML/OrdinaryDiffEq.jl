using OrdinaryDiffEq, ForwardDiff, Test

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/4388
# u' = 1, u(0) = 0 fires `u - q = 0` at T = q, so dT/dq = 1. The discrete event time has
# dT/dq = 1/û'(T) for the dense output û, whose slope error is ~1e-12 for Vern9 at these tolerances.
function event_time(q, alg; iip = false)
    cb = ContinuousCallback((u, t, integ) -> (iip ? u[1] : u) - q, terminate!)
    tspan = (zero(q), 10 * one(q))
    prob = iip ?
        ODEProblem((du, u, p, t) -> (du[1] = one(eltype(du))), [zero(q)], tspan) :
        ODEProblem((u, p, t) -> one(typeof(q)), zero(q), tspan)
    sol = solve(
        prob, alg; reltol = 1.0e-14, abstol = 1.0e-14,
        callback = cb, save_everystep = false
    )
    return sol.t[end]
end

# u' = p u, u(0) = 1 fires `u - 2 = 0` at T = log(2)/p, so dT/dp = -log(2)/p^2.
function growth_event_time(p, alg)
    cb = ContinuousCallback((u, t, integ) -> u - 2, terminate!)
    prob = ODEProblem((u, p, t) -> p * u, one(p), (zero(p), 10 * one(p)), p)
    sol = solve(prob, alg; reltol = 1.0e-13, abstol = 1.0e-13, callback = cb)
    return sol.t[end]
end

@testset "ForwardDiff event time through ContinuousCallback" begin
    for alg in (Tsit5(), Vern9()), q in (3.0, 1.7, 5.3), iip in (false, true)
        @test event_time(q, alg; iip) ≈ q rtol = 1.0e-12
        @test ForwardDiff.derivative(q -> event_time(q, alg; iip), q) ≈ 1 rtol = 1.0e-10
    end
    for alg in (Tsit5(), Vern9()), p in (0.7, 1.3)
        @test ForwardDiff.derivative(p -> growth_event_time(p, alg), p) ≈ -log(2) / p^2 rtol = 1.0e-9
    end
end
