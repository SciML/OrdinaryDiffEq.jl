@inline function initialize!(integrator,cache::SSPRK22ConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

@inline function perform_step!(integrator,cache::SSPRK22ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  k = integrator.fsalfirst
  tmp = uprev + dt*k
  k = f(t+dt,tmp)
  u = (uprev + tmp + dt*k) / 2
  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::SSPRK22Cache,f=integrator.f)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@inline function perform_step!(integrator,cache::SSPRK22Cache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k,du,tmp,fsalfirst = cache
  for i in uidx
    tmp[i] = @muladd uprev[i] + dt*integrator.fsalfirst[i]
  end
  f(t+dt,tmp,k)
  for i in uidx
    u[i] = (uprev[i] + tmp[i] + dt*k[i]) / 2
  end
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end


@inline function initialize!(integrator,cache::SSPRK33ConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

@inline function perform_step!(integrator,cache::SSPRK33ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  k = integrator.fsalfirst
  tmp = uprev + dt*k
  k = f(t+dt,tmp)
  tmp = (3*uprev + tmp + dt*k) / 4
  k = f(t+dt/2,tmp)
  u = (uprev + 2*tmp + 2*dt*k) / 3
  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::SSPRK33Cache,f=integrator.f)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@inline function perform_step!(integrator,cache::SSPRK33Cache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k,du,tmp,fsalfirst = cache
  for i in uidx
    tmp[i] = @muladd uprev[i] + dt*integrator.fsalfirst[i]
  end
  f(t+dt,tmp,k)
  for i in uidx
    tmp[i] = (3*uprev[i] + tmp[i] + dt*k[i]) / 4
  end
  f(t+dt/2,tmp,k)
  for i in uidx
    u[i] = (uprev[i] + 2*tmp[i] + 2*dt*k[i]) / 3
  end
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end
