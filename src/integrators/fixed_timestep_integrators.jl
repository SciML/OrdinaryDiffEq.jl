@inline function initialize!(integrator,cache::EulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
end

@inline function perform_step!(integrator::ODEIntegrator,cache::EulerConstantCache)
  @unpack t,dt,uprev,u,f,k = integrator
  k = integrator.fsalfirst
  u = muladd(dt,k,uprev)
  k = f(t+dt,u) # For the interpolation, needs k at the updated point
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::EulerCache)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@inline function perform_step!(integrator::ODEIntegrator,cache::EulerCache)
  @unpack t,dt,uprev,u,f,k = integrator
  uidx = eachindex(integrator.uprev)
  for i in uidx
    u[i] = muladd(dt,integrator.fsalfirst[i],uprev[i])
  end
  f(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
  @pack integrator = t,dt,u,k
end

@inline function initialize!(integrator,cache::MidpointConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

@inline function perform_step!(integrator::ODEIntegrator,cache::MidpointConstantCache)
  @unpack t,dt,uprev,u,f,k = integrator
  halfdt = dt/2
  k = integrator.fsalfirst
  k = f(t+halfdt,uprev+halfdt*k)
  u = uprev + dt*k
  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::MidpointCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@inline function perform_step!(integrator::ODEIntegrator,cache::MidpointCache)
  @unpack t,dt,uprev,u,f,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k,du,utilde,fsalfirst = cache
  halfdt = dt/2
  for i in uidx
    utilde[i] = muladd(halfdt,integrator.fsalfirst[i],uprev[i])
  end
  f(t+halfdt,utilde,du)
  for i in uidx
    u[i] = muladd(dt,du[i],uprev[i])
  end
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::RK4ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

@inline function perform_step!(integrator::ODEIntegrator,cache::RK4ConstantCache)
  @unpack t,dt,uprev,u,f,k = integrator
  halfdt = dt/2
  k₁ =integrator.fsalfirst
  ttmp = t+halfdt
  k₂ = f(ttmp,muladd(halfdt,k₁,uprev))
  k₃ = f(ttmp,muladd(halfdt,k₂,uprev))
  k₄ = f(t+dt,muladd(dt,k₃,uprev))
  u = muladd(dt/6,muladd(2,(k₂ + k₃),k₁+k₄),uprev)
  k = f(t+dt,u)
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::RK4Cache)
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # pre-start FSAL
end

@inline function perform_step!(integrator::ODEIntegrator,cache::RK4Cache)
  @unpack t,dt,uprev,u,f,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  k₁ = fsalfirst
  halfdt = dt/2
  ttmp = t+halfdt
  for i in uidx
    tmp[i] = muladd(halfdt,k₁[i],uprev[i])
  end
  f(ttmp,tmp,k₂)
  for i in uidx
    tmp[i] = muladd(halfdt,k₂[i],uprev[i])
  end
  f(ttmp,tmp,k₃)
  for i in uidx
    tmp[i] = muladd(dt,k₃[i],uprev[i])
  end
  f(t+dt,tmp,k₄)
  for i in uidx
    u[i] = muladd(dt/6,muladd(2,(k₂[i] + k₃[i]),k₁[i] + k₄[i]),uprev[i])
  end
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end
