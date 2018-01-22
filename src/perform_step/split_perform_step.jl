function initialize!(integrator,cache::SplitEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f.f1(integrator.t,integrator.uprev) + integrator.f.f2(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::SplitEulerConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  u = @. uprev + dt*integrator.fsalfirst
  integrator.fsallast = f.f1(t+dt,u) + f.f2(t+dt,u)  # For the interpolation, needs k at the updated point
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::SplitEulerCache)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f.f1(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
  integrator.f.f2(integrator.t,integrator.uprev,cache.tmp) # For the interpolation, needs k at the updated point
  integrator.fsalfirst .+= cache.tmp
end

@muladd function perform_step!(integrator,cache::SplitEulerCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @. u = uprev + dt*integrator.fsalfirst
  f.f1(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
  f.f2(t+dt,u,cache.tmp) # For the interpolation, needs k at the updated point
  integrator.fsallast .+= cache.tmp
end
