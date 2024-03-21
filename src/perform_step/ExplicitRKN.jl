# three steps:
# 1) cache
#   a)struct that contains all of the type needed for the method
#   b) alg_cache which caches the algorithm
# 2) intialize!: Set the values of the integrator using the cache.
# 3) perform_step: perform one iteration of the nystrom method
#       set the nodes to take the slopes at
#       derive each k value
#       find the new u and u' values.

using OrdinaryDiffEq
import OrdinaryDiffEq:
                       OrdinaryDiffEqAlgorithm, OrdinaryDiffEqMutableCache,
                       OrdinaryDiffEqConstantCache,
                       alg_order, alg_cache, initialize!, perform_step!, trivial_limiter!,
                       constvalue,
                       @muladd, @unpack, @cache, @..

struct RKN4_ALG{StageLimiter, StepLimiter} <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
end
RKN4_ALG(stage_limiter! = trivial_limiter!) = RKN4_ALG(stage_limiter!, trivial_limiter!)
export RKN_ALG
alg_order(alg::RKN4_ALG) = 3


@cache struct RKN4Cache{uType, rateType, reducedRateType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k₂::reducedRateType
    k₃::reducedRateType
    k::rateType
    tmp::uType
end

function alg_cache(alg::RKN4, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    k₁ = zero(rate_prototype)
    k₂ = zero(reduced_rate_prototype)
    k₃ = zero(reduced_rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    ExplicitRKN4Cache(u, uprev, k₁, k₂, k₃, k, tmp)
end

struct RKNconstantCache <: OrdinaryDiffEqConstantCache end

function initialize!(integrator, cache::RKN4Cache)
    @unpack fsalfirst, k = cache
    duprev, uprev = integrator.uprev.x
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.k[1].x[1], duprev, uprev, integrator.p, integrator.t)
    integrator.f.f2(integrator.k[1].x[2], duprev, uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
end

@muladd function perform_step!(integrator, cache::RKNconstantCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    
    #define dt values
    halfdt = dt/2
    dtsq = dt^2
    eightdtsq = dtsq/8
    halfdtsq = dtsq/2
    sixthdtsq = dtsq/6
    sixthdt = dt/6
    ttmp = t + halfdt

    #perform operations to find k values
    k₁ = integrator.fsalfirst.x[1]
    ku = uprev + halfdt * duprev + eightdtsq * k₁
    kdu = duprev + halfdt * k₁

    k₂ = f.f1(kdu, ku, p, ttmp)
    ku = uprev + dt * duprev + halfdtsq * k₂
    kdu = duprev + dt * k₂

    k₃ = f.f1(kdu, ku, p, t + dt)
    ku = uprev + dt * duprev + eightdtsq * k₃
    kdu = duprev + dt * k₃

    #perform final calculations to determine new y and y'.
    u = uprev + sixthdtsq* (1*k₁ + 2*k₂ + 0*k₃) + dt * duprev
    du = duprev + sixthdt * (1*k₁ + 4*k₂ + 1*k₃)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 3
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

f = ODEFunction((u, p, t) -> 1.01u,
    analytic = (u0, p, t) -> u0 * exp(1.01t))
prob = ODEProblem(f, 1.01, (0.0, 1.0))
sol = solve(prob, RKN4_ALG(), dt = 0.1)

using Plots
plot(sol)
plot(sol, denseplot = false, plot_analytic = true)

using DiffEqDevTools
dts = (1 / 2) .^ (8:-1:1)
sim = test_convergence(dts, prob, RKN4_ALG())
sim.𝒪est[:final]
plot(sim)

# Example of a good one!
sim = test_convergence(dts, prob, BS3())
sim.𝒪est[:final]
plot(sim)

#iip test
f = ODEFunction((du, u, p, t) -> (du .= 1.01 .* u),
    analytic = (u0, p, t) -> u0 * exp(1.01t))
prob = ODEProblem(f, [1.01], (0.0, 1.0))
sol = solve(prob, RKN4_ALG(), dt = 0.1)

plot(sol)
plot(sol, denseplot = false, plot_analytic = true)

dts = (1 / 2) .^ (8:-1:1)
sim = test_convergence(dts, prob, RKN4_ALG())
sim.𝒪est[:final]
plot(sim)

# Example of a good one!
sim = test_convergence(dts, prob, BS3())
sim.𝒪est[:final]
plot(sim)


