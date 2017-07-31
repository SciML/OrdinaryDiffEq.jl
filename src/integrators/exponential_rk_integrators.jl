@inline function initialize!(integrator,cache::LawsonEulerConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  rtmp = f[2]
  integrator.fsalfirst = rtmp # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = zero(integrator.fsalfirst)
end

@inline @muladd function perform_step!(integrator,cache::LawsonEulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  rtmp = integrator.fsalfirst
  A = f[1]
  u = expm(dt*A)*(@. uprev + dt*rtmp)
  rtmp = f[2](t+dt,u)
  k = A*u .+ rtmp # For the interpolation, needs k at the updated point
  integrator.fsallast = rtmp
  integrator.k[1] = integrator.fsalfirst # this is wrong, since it's just rtmp. Should fsal this value though
  integrator.k[2] = k
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::LawsonEulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst,rtmp = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = rtmp
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = fsalfirst # this is wrong, since it's just rtmp. Should fsal this value though
  integrator.k[2] = k
  A = f[1]
  A_mul_B!(cache.k,A,integrator.u)
  f[2](integrator.t,integrator.uprev,rtmp) # For the interpolation, needs k at the updated point
  @. integrator.fsalfirst = cache.k + rtmp
end

@inline @muladd function perform_step!(integrator,cache::LawsonEulerCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack k,rtmp,tmp = cache
  A = f[1]
  M = expm(dt*A)
  @. tmp = uprev + dt*integrator.fsalfirst
  A_mul_B!(u,M,tmp)
  A_mul_B!(tmp,A,u)
  f[2](t+dt,u,rtmp)
  @. k = tmp +  rtmp
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::NorsettEulerConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  rtmp = f[2](integrator.t,integrator.uprev)
  integrator.fsalfirst = rtmp # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = zero(integrator.fsalfirst)
end

@inline function perform_step!(integrator,cache::NorsettEulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  rtmp = integrator.fsalfirst
  A = f[1]
  u = uprev .+ ((expm(dt*A)-I)/A)*(A*uprev .+ rtmp)
  rtmp = f[2](t+dt,u)
  k = A*u .+ rtmp # For the interpolation, needs k at the updated point
  integrator.fsallast = rtmp
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = k
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::NorsettEulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst,rtmp = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = rtmp
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = fsalfirst
  integrator.k[2] = k
  A = f[1](integrator.t,integrator.u,rtmp)
  A_mul_B!(cache.k,A,integrator.u)
  f[2](integrator.t,integrator.uprev,rtmp) # For the interpolation, needs k at the updated point
  @. integrator.fsalfirst = cache.k + rtmp
end

@inline function perform_step!(integrator,cache::NorsettEulerCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack k,rtmp,tmp = cache
  A = f[1]
  M = ((expm(dt*A)-I)/A)
  A_mul_B!(tmp,A,uprev)
  tmp .+= rtmp
  A_mul_B!(rtmp,M,tmp)
  @. u = uprev + rtmp
  A_mul_B!(tmp,A,u)
  f[2](t+dt,u,rtmp)
  @. k = tmp +  rtmp
  @pack integrator = t,dt,u
end
