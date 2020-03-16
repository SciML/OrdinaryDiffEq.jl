using OrdinaryDiffEq
import OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqConstantCache, OrdinaryDiffEqMutableCache,
      alg_order, alg_cache, initialize!, perform_step!, @muladd, @unpack, @pack!,
      constvalue

struct LDDRK <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end
export LDDRK
alg_order(alg::LDDRK) = 4

struct LDDRK_Cache <: OrdinaryDiffEqConstantCache
   α1
   α2
   α3
   α4
   α5
   β1
   β2
   β3
   β4
   β5
   c1
   c2
   c3
   c4
   c5
end

 function LDDRK_Cache(T)
   α1 = T(0.0)
   α2 = T(-0.6913065)
   α3 = T(-2.655155)
   α4 = T(-0.8147688)
   α5 = T(-0.66865870)
   β1 = T(0.1)
   β2 = T(0.75)
   β3 = T(0.7)
   β4 = T(0.479313)
   β5 = T(0.310392)
   c1 = T(0.0)
   c2 = T(0.1)
   c3 = T(0.3315201)
   c4 = T(0.4577796)
   c5 = T(0.86665284)
   LDDRK_Cache(α1, α2, α3, α4, α5, β1, β2, β3, β4, β5, c1, c2, c3, c4, c5)
end

alg_cache(alg::LDDRK,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) = LDDRK_Cache(tTypeNoUnits)

function initialize!(integrator, cache::LDDRK_Cache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast

end

@muladd function perform_step!(integrator, cache::LDDRK_Cache, repeat_step=false)
  @unpack t,dt,uprev,u, fsallast, f,p = integrator
  @unpack α1, α2, α3, α4, α5, β1, β2, β3, β4, β5, c1, c2, c3, c4, c5 = cache

  γ = fsallast

  γ = γ*α1 + dt*f(u, p, t + c1*dt)
  u = u + β1*γ

  γ = γ*α2 + dt*f(u, p, t + c2*dt)
  u = u + β2*γ

  γ = γ*α3 + dt*f(u, p, t + c3*dt)
  u = u + β3*γ

  γ = γ*α4 + dt*f(u, p, t + c4*dt)
  u = u + β4*γ

  γ = γ*α5 + dt*f(u, p, t + c5*dt)
  u = u + β5*γ

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

f = ODEFunction((u,p,t)->-u,
            analytic = (u0,p,t) -> u0*exp(-t))
prob = ODEProblem(f,1.01,(0.0,1.0))

using Plots
using DiffEqDevTools
dts = (1/2) .^ (8:-1:1)
sim = test_convergence(dts,prob,LDDRK())
plot(sim)
