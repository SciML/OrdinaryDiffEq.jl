function initialize!(integrator,cache::AB3ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::AB3ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  k1 = integrator.fsalfirst
  cnt = length(integrator.sol.t)
  if cnt == 1 || cnt == 2
    ttmp = t + (2/3)*dt
    ralk2 = f(uprev + (2/3)*dt*k1, p, ttmp)
    u = uprev + (dt/4)*(k1 + 3*ralk2)
  else
    u1 = integrator.sol(t-dt)
    u2 = integrator.sol(t-2*dt)
    k2 = f(u1, p, t-dt)
    k3 = f(u2, p, t-2*dt)
    u  = uprev + (dt/12)*(23*k1 - 16*k2 + 5*k3)
  end
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::AB3Cache)
  @unpack tmp,fsalfirst,k2,k3,k = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # pre-start FSAL
end

@muladd function perform_step!(integrator,cache::AB3Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,fsalfirst,k2,k3,ralk2,u1,u2,k = cache
  k1 = integrator.fsalfirst
  cnt = length(integrator.sol.t)
  if cnt == 1 || cnt == 2
    @. ttmp = t + (2/3)*dt
    f(ralk2, uprev + (2/3)*dt*k1, p, ttmp)
    @. u = uprev + (dt/4)*(k1 + 3*k2)
  else
    u1 = integrator.sol(t-dt)
    u2 = integrator.sol(t-2*dt)
    f(k2, u1, p, t-dt)
    f(k3, u2, p, t-2*dt)
    @. u  = uprev + (dt/12)*(23*k1 - 16*k2 + 5*k3)
  end
  f(k, u, p, t+dt)
end

function initialize!(integrator,cache::ABM32ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::ABM32ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  k1 = integrator.fsalfirst
  cnt = length(integrator.sol.t)
  if cnt == 1
    ttmp = t + (2/3)*dt
    ralk2 = f(uprev + (2/3)*dt*k1, p, ttmp)
    u = uprev + (dt/4)*(k1 + 3*ralk2)
  else
    perform_step!(integrator,AB3ConstantCache())
    u2 = integrator.sol(t-dt)
    k2 = integrator.fsallast
    k3 = f(u2, p, t-dt)
    u = uprev + (dt/12)*(5*k2 + 8*k1 - k3)
  end
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::ABM32Cache)
  @unpack tmp,fsalfirst,k2,k3,k = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # pre-start FSAL
end

@muladd function perform_step!(integrator,cache::ABM32Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,fsalfirst,k2,k3,ralk2,u2,k = cache
  k1 = integrator.fsalfirst
  cnt = length(integrator.sol.t)
  if cnt == 1
    @. ttmp = t + (2/3)*dt
    f(ralk2, uprev + (2/3)*dt*k1, p, ttmp)
    @. u = uprev + (dt/4)*(k1 + 3*ralk2)
  else
    perform_step!(integrator,integrator.cache.onestep_cache)
    u2 = integrator.sol(t-dt)
    k2 = integrator.fsallast
    f(k3, u2, p, t-dt)
    @. u = uprev + (dt/12)*(5*k2 + 8*k1 - k3)
  end
  f(k, u, p, t+dt)
end