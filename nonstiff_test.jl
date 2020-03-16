using OrdinaryDiffEq
import OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqConstantCache, OrdinaryDiffEqMutableCache,
      alg_order, alg_cache, initialize!, perform_step!, @muladd, @unpack, @pack!,
      constvalue

struct LDDRK <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end
export LDDRK
alg_order(alg::LDDRK) = 4

mutable struct LDDRKcache <: OrdinaryDiffEqMutableCache
   alpha1
   alpha2
   alpha3
   alpha4
   alpha5
   beta1
   beta2
   beta3
   beta4
   beta5
   c1
   c2
   c3
   c4
   c5
   uCache
   wCache
end

 function LDDRKcache(T)
   alpha1 = T(0.0)
   alpha2 = T(-0.6913065)
   alpha3 = T(-2.655155)
   alpha4 = T(-0.8147688)
   alpha5 = T(-0.66865870)
   beta1 = T(0.1)
   beta2 = T(0.75)
   beta3 = T(0.7)
   beta4 = T(0.479313)
   beta5 = T(0.310392)
   c1 = T(0.0)
   c2 = T(0.1)
   c3 = T(0.3315201)
   c4 = T(0.4577796)
   c5 = T(0.86665284)
   uCache = T(0)
   wCache = T(0)
   LDDRKcache(alpha1, alpha2, alpha3, alpha4, alpha5, beta1, beta2, beta3, beta4, beta5, c1, c2, c3, c4, c5, uCache, wCache)
end

alg_cache(alg::LDDRK,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) = LDDRKcache(tTypeNoUnits)

function initialize!(integrator, cache::LDDRKcache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uCache = integrator.uprev
  wCache = integrator.f(integrator.uprev, integrator.p, integrator.t)
  @pack! cache = uCache, wCache

end

@muladd function perform_step!(integrator, cache::LDDRKcache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack alpha1, alpha2, alpha3, alpha4, alpha5, beta1, beta2, beta3, beta4, beta5, c1, c2, c3, c4, c5, uCache, wCache = cache
  alpha = [alpha1, alpha2, alpha3, alpha4, alpha5]
  beta = [beta1, beta2, beta3, beta4, beta5]
  c = [c1, c2, c3, c4, c5]
  for i in range(1,5)
    wCache = wCache*alpha[i] + dt*f(uCache, p, t + c[i]*dt)
    uCache = uCache + beta[i]*wCache
  end
  @pack! cache = uCache, wCache
  u = uCache
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end
