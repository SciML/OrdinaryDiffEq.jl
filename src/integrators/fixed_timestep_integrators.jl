@inline function initialize!(integrator,cache::DiscreteConstantCache,f=integrator.f)
  integrator.kshortsize = 0
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

@inline function initialize!(integrator,cache::DiscreteCache,f=integrator.f)
  integrator.kshortsize = 0
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

@inline function perform_step!(integrator,cache::DiscreteConstantCache,f=integrator.f)
  if discrete_apply_map(integrator.alg)
    if discrete_scale_by_time(integrator.alg)
      integrator.u = integrator.uprev .+ integrator.dt.*f(integrator.t+integrator.dt,integrator.uprev)
    else
      integrator.u = f(integrator.t+integrator.dt,integrator.uprev)
    end
  end
end

@inline function perform_step!(integrator,cache::DiscreteCache,f=integrator.f)
  @unpack u,uprev,dt,t = integrator
  @unpack du = cache
  if discrete_apply_map(integrator.alg)
    if discrete_scale_by_time(integrator.alg)
      f(t+dt,uprev,du)
      @. u = @muladd uprev + dt*du
    else
      f(t+dt,uprev,u)
    end
    if typeof(uprev) <: DEDataArray # Needs to get the fields, since updated uprev
      copy_non_array_fields!(u,uprev)
    end
  end
end

@inline function initialize!(integrator,cache::EulerConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
end

@inline function perform_step!(integrator,cache::EulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  k = integrator.fsalfirst
  u = muladd.(dt,k,uprev)
  k = f(t+dt,u) # For the interpolation, needs k at the updated point
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::EulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@inline function perform_step!(integrator,cache::EulerCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @. u = muladd(dt,integrator.fsalfirst,uprev)
  f(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
  @pack integrator = t,dt,u,k
end

@inline function initialize!(integrator,cache::MidpointConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

@inline function perform_step!(integrator,cache::MidpointConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  halfdt = dt/2
  k = integrator.fsalfirst
  k = f(t+halfdt,uprev.+halfdt.*k)
  u = uprev .+ dt.*k
  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::MidpointCache,f=integrator.f)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@inline function perform_step!(integrator,cache::MidpointCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k,du,tmp,fsalfirst = cache
  halfdt = dt/2
  @. tmp = muladd(halfdt,integrator.fsalfirst,uprev)
  f(t+halfdt,tmp,du)
  @. u = muladd(dt,du,uprev)
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::RK4ConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

@inline function perform_step!(integrator,cache::RK4ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  halfdt = dt/2
  k₁ =integrator.fsalfirst
  ttmp = t+halfdt
  k₂ = f(ttmp,muladd.(halfdt,k₁,uprev))
  k₃ = f(ttmp,muladd.(halfdt,k₂,uprev))
  k₄ = f(t+dt,muladd.(dt,k₃,uprev))
  u = muladd.(dt/6,muladd.(2,(k₂ .+ k₃),k₁.+k₄),uprev)
  k = f(t+dt,u)
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::RK4Cache,f=integrator.f)
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # pre-start FSAL
end

@inline function perform_step!(integrator,cache::RK4Cache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  k₁ = fsalfirst
  halfdt = dt/2
  ttmp = t+halfdt
  @. tmp = muladd(halfdt,k₁,uprev)
  f(ttmp,tmp,k₂)
  @. tmp = muladd(halfdt,k₂,uprev)
  f(ttmp,tmp,k₃)
  @. tmp = muladd(dt,k₃,uprev)
  f(t+dt,tmp,k₄)
  @. u = muladd(dt/6,muladd(2,(k₂ + k₃),k₁ + k₄),uprev)
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end
