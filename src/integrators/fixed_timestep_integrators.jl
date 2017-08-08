function initialize!(integrator,cache::DiscreteConstantCache,f=integrator.f)
  integrator.kshortsize = 0
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

function perform_step!(integrator,cache::DiscreteConstantCache,f=integrator.f)
  if discrete_apply_map(integrator.alg)
    if discrete_scale_by_time(integrator.alg)
      @muladd integrator.u = integrator.uprev .+ integrator.dt.*f(integrator.t+integrator.dt,integrator.uprev)
    else
      integrator.u = f(integrator.t+integrator.dt,integrator.uprev)
    end
  end
end

function initialize!(integrator,cache::DiscreteCache,f=integrator.f)
  integrator.kshortsize = 0
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

function perform_step!(integrator,cache::DiscreteCache,f=integrator.f)
  @unpack u,uprev,dt,t = integrator
  @unpack du = cache
  if discrete_apply_map(integrator.alg)
    if discrete_scale_by_time(integrator.alg)
      f(t+dt,uprev,du)
      @muladd @. u = uprev + dt*du
    else
      f(t+dt,uprev,u)
    end
    if typeof(uprev) <: DEDataArray # Needs to get the fields, since updated uprev
      copy_fields!(u,uprev)
    end
  end
end

function initialize!(integrator,cache::EulerConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::EulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @muladd u = @. uprev + dt*integrator.fsalfirst
  k = f(t+dt,u) # For the interpolation, needs k at the updated point
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::EulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

function perform_step!(integrator,cache::EulerCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @muladd @. u = uprev + dt*integrator.fsalfirst
  f(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
end

function initialize!(integrator,cache::Union{HeunConstantCache,RalstonConstantCache},f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::Union{HeunConstantCache,RalstonConstantCache},f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  fsalfirst = integrator.fsalfirst

  if typeof(cache) <: HeunConstantCache
      a = dt
  else # Ralston
      a = 3dt/4
  end

  @muladd tmp = @. uprev + a*fsalfirst

  k2 = f(t+a,tmp)

  if typeof(cache) <: HeunConstantCache
      @muladd u = @. uprev + (dt/2)*(fsalfirst + k2)
  else
      @muladd u = @. uprev + (dt/3)*fsalfirst + (2dt/3)*k2
  end

  if integrator.opts.adaptive
      if typeof(cache) <: HeunConstantCache
          @muladd utilde = @. (dt/2)*(k2 - fsalfirst)
      else
          @muladd utilde = @. (2dt/3)*(k2 - fsalfirst)
      end

      tmp = @. utilde/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
      integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  k = f(t+dt,u)
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::Union{HeunCache,RalstonCache},f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

function perform_step!(integrator,cache::Union{HeunCache,RalstonCache},f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack fsalfirst,k,utilde = cache

  if typeof(cache) <: HeunCache
      a = dt
  else # Ralston
      a = 3dt/4
  end

  @muladd @. utilde = uprev + a*fsalfirst

  f(t+a,utilde,k)

  if typeof(cache) <: HeunCache
      @muladd @. u = uprev + (dt/2)*(fsalfirst + k)
  else
      @muladd @. u = uprev + (dt/3)*fsalfirst + (2dt/3)*k
  end

  if integrator.opts.adaptive
      if typeof(cache) <: HeunCache
          @muladd @. utilde = (dt/2)*(k - fsalfirst)
      else
          @muladd @. utilde = (2dt/3)*(k - fsalfirst)
      end

      @. utilde = utilde/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
      integrator.EEst = integrator.opts.internalnorm(utilde)
  end
  f(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
end

function initialize!(integrator,cache::MidpointConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::MidpointConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  halfdt = dt/2
  k = f(t+halfdt, @. uprev + halfdt*integrator.fsalfirst)
  u = @. uprev + dt*k
  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  if integrator.opts.adaptive
      utilde = @. dt*(integrator.fsalfirst - k)
      tmp = @. utilde/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
      integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::MidpointCache,f=integrator.f)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::MidpointCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack tmp,k,fsalfirst,utilde = cache
  halfdt = dt/2
  @. tmp = uprev + halfdt*fsalfirst
  f(t+halfdt,tmp,k)
  @. u = uprev + dt*k
  if integrator.opts.adaptive
      @. utilde = dt*(fsalfirst - k)
      @. utilde = (utilde)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
      integrator.EEst = integrator.opts.internalnorm(utilde)
  end
  f(t+dt,u,k)
end

function initialize!(integrator,cache::RK4ConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::RK4ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  halfdt = dt/2
  k₁ =integrator.fsalfirst
  ttmp = t+halfdt
  k₂ = f(ttmp, @. uprev + halfdt*k₁)
  k₃ = f(ttmp, @. uprev + halfdt*k₂)
  k₄ = f(t+dt, @. uprev + dt*k₃)
  u = @. uprev + (dt/6)*(2*(k₂ + k₃) + (k₁+k₄))
  integrator.fsallast = f(t+dt,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::RK4Cache,f=integrator.f)
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # pre-start FSAL
end

@muladd function perform_step!(integrator,cache::RK4Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  k₁ = fsalfirst
  halfdt = dt/2
  ttmp = t+halfdt
  @. tmp = uprev + halfdt*k₁
  f(ttmp,tmp,k₂)
  @. tmp = uprev + halfdt*k₂
  f(ttmp,tmp,k₃)
  @. tmp = uprev + dt*k₃
  f(t+dt,tmp,k₄)
  @. u = uprev + (dt/6)*(2*(k₂ + k₃) + (k₁ + k₄))
  f(t+dt,u,k)
end
