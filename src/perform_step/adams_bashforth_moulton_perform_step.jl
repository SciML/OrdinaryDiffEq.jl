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
  @unpack k2, k3 = cache
  k1 = integrator.fsalfirst
  cnt = integrator.iter
  if cnt == 1 || cnt == 2
    ttmp = t + (2/3)*dt
    ralk2 = f(uprev + (2/3)*dt*k1, p, ttmp)       #Ralston Method
    u = uprev + (dt/4)*(k1 + 3*ralk2)
    if cnt == 1
      k3 = k1
    else
      k2 = k1
    end
  else
    u  = uprev + (dt/12)*(23*k1 - 16*k2 + 5*k3)
    k3 = k2
    k2 = k1
  end
  cache.k2 = k2
  cache.k3 = k3
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::AB3Cache)
  @unpack tmp,fsalfirst,k = cache
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
  @unpack tmp,fsalfirst,k2,k3,ralk2,k = cache
  k1 = integrator.fsalfirst
  cnt = integrator.iter
  if cnt == 1 || cnt == 2
    ttmp = t + (2/3)*dt
    @. tmp = uprev + (2/3)*dt*k1
    f(ralk2, tmp, p, ttmp)    
    @. u = uprev + (dt/4)*(k1 + 3*ralk2)        #Ralston Method
    if cnt == 1
      cache.k3 = copy(k1)
    else
      cache.k2 = copy(k1)
    end
  else
    @. u  = uprev + (dt/12)*(23*k1 - 16*k2 + 5*k3)
    cache.k3 = copy(k2)
    cache.k2 = copy(k1)
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
  @unpack k2,k3 = cache
  k1 = integrator.fsalfirst
  cnt = integrator.iter
  if cnt == 1
    ttmp = t + (2/3)*dt
    ralk2 = f(uprev + (2/3)*dt*k1, p, ttmp)     #Ralston Method
    u = uprev + (dt/4)*(k1 + 3*ralk2)
    k2 = k1
  else
    perform_step!(integrator, AB3ConstantCache(k2,k3))
    k = integrator.fsallast
    u = uprev + (dt/12)*(5*k + 8*k1 - k2)
    k3 = k2
    k2 = k1
  end
  cache.k2 = k2
  cache.k3 = k3
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::ABM32Cache)
  @unpack fsalfirst,k = cache
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
  @unpack tmp,fsalfirst,k2,k3,ralk2,k = cache
  k1 = integrator.fsalfirst
  cnt = integrator.iter
  if cnt == 1
    ttmp = t + (2/3)*dt
    @. tmp = uprev + (2/3)*dt*k1
    f(ralk2, tmp, p, ttmp)
    @. u = uprev + (dt/4)*(k1 + 3*ralk2)       #Ralston Method
    cache.k2 = copy(k1)
  else
    perform_step!(integrator, AB3Cache(u,uprev,fsalfirst,k2,k3,ralk2,k,tmp))
    k = integrator.fsallast
    @. u = uprev + (dt/12)*(5*k + 8*k1 - k2)
    cache.k3 = copy(k2)
    cache.k2 = copy(k1)
  end
  f(k, u, p, t+dt)
end
