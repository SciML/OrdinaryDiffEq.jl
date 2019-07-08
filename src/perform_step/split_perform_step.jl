function initialize!(integrator,cache::SplitEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst =
  integrator.f.f1(integrator.uprev, integrator.t, integrator) +
  integrator.f.f2(integrator.uprev, integrator.t, integrator) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.destats.nf2 += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::SplitEulerConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  u = @.. uprev + dt*integrator.fsalfirst
  integrator.fsallast = f.f1(u, t+dt, integrator) + f.f2(u, t+dt, integrator)  # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
  integrator.destats.nf2 += 1
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
  integrator.f.f1(integrator.fsalfirst, integrator.uprev, integrator.t, integrator) # For the interpolation, needs k at the updated point
  integrator.f.f2(cache.tmp, integrator.uprev, integrator.t, integrator) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
  integrator.destats.nf2 += 1
  integrator.fsalfirst .+= cache.tmp
end

@muladd function perform_step!(integrator,cache::SplitEulerCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @.. u = uprev + dt*integrator.fsalfirst
  f.f1(integrator.fsallast, u, t+dt, integrator) # For the interpolation, needs k at the updated point
  f.f2(cache.tmp, u, t+dt, integrator) # For the interpolation, needs k at the updated point
  integrator.destats.nf2 += 1
  integrator.destats.nf += 1
  integrator.fsallast .+= cache.tmp
end
