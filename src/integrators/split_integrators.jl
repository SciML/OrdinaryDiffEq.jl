@inline function initialize!(integrator,cache::SplitEulerConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f[1](integrator.t,integrator.uprev) + f[2](integrator.t,integrator.uprev) # Pre-start fsal
end

@inline function perform_step!(integrator,cache::SplitEulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  k = integrator.fsalfirst
  u = muladd(dt,k,uprev)
  k = f[1](t+dt,u) + f[2](t+dt,u)  # For the interpolation, needs k at the updated point
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::SplitEulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f[1](integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
  f[2](integrator.t,integrator.uprev,integrator.cache.tmp) # For the interpolation, needs k at the updated point
  for i in eachindex(integrator.u)
    integrator.fsalfirst[i] += cache.tmp[i]
  end
end

@inline function perform_step!(integrator,cache::SplitEulerCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  for i in uidx
    u[i] = muladd(dt,integrator.fsalfirst[i],uprev[i])
  end
  f[1](t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
  f[2](t+dt,u,integrator.cache.tmp) # For the interpolation, needs k at the updated point
  for i in uidx
    integrator.fsallast[i] += cache.tmp[i]
  end
  @pack integrator = t,dt,u,k
end
