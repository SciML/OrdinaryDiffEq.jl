@inline function initialize!(integrator,cache::SplitEulerConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f[1](integrator.t,integrator.uprev) + f[2](integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline @muladd function perform_step!(integrator,cache::SplitEulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  u = @. uprev + dt*integrator.fsalfirst
  integrator.fsallast = f[1](t+dt,u) + f[2](t+dt,u)  # For the interpolation, needs k at the updated point
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
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
  f[2](integrator.t,integrator.uprev,cache.tmp) # For the interpolation, needs k at the updated point
  integrator.fsalfirst .+= cache.tmp
end

@inline @muladd function perform_step!(integrator,cache::SplitEulerCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @. u = uprev + dt*integrator.fsalfirst
  f[1](t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
  f[2](t+dt,u,cache.tmp) # For the interpolation, needs k at the updated point
  integrator.fsallast .+= cache.tmp
end
