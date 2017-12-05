function initialize!(integrator,cache::DiscreteConstantCache)
  integrator.kshortsize = 0
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
end

function perform_step!(integrator,cache::DiscreteConstantCache,repeat_step=false)
  if discrete_apply_map(integrator.alg)
    if discrete_scale_by_time(integrator.alg)
      @muladd integrator.u = integrator.uprev .+ integrator.dt.*integrator.f(integrator.t+integrator.dt,integrator.uprev)
    else
      integrator.u = integrator.f(integrator.t+integrator.dt,integrator.uprev)
    end
  end
end

function initialize!(integrator,cache::DiscreteCache)
  integrator.kshortsize = 0
  resize!(integrator.k, integrator.kshortsize)
end

function perform_step!(integrator,cache::DiscreteCache,repeat_step=false)
  @unpack u,uprev,dt,t,f = integrator
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

function initialize!(integrator,cache::EulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::EulerConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @muladd u = @. uprev + dt*integrator.fsalfirst
  k = f(t+dt,u) # For the interpolation, needs k at the updated point
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::EulerCache)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

function perform_step!(integrator,cache::EulerCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @muladd @. u = uprev + dt*integrator.fsalfirst
  f(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
end

function initialize!(integrator,cache::Union{HeunConstantCache,RalstonConstantCache})
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::Union{HeunConstantCache,RalstonConstantCache},repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
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

      tmp = @. utilde/(integrator.opts.abstol+max(integrator.opts.internalnorm(uprev),
                        integrator.opts.internalnorm(u))*integrator.opts.reltol)
      integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  k = f(t+dt,u)
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::Union{HeunCache,RalstonCache})
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

function perform_step!(integrator,cache::Union{HeunCache,RalstonCache},repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
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

      @. utilde = utilde/(integrator.opts.abstol+max(integrator.opts.internalnorm(uprev),
                          integrator.opts.internalnorm(u))*integrator.opts.reltol)
      integrator.EEst = integrator.opts.internalnorm(utilde)
  end
  f(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
end

function initialize!(integrator,cache::MidpointConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::MidpointConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  halfdt = dt/2
  k = f(t+halfdt, @. uprev + halfdt*integrator.fsalfirst)
  u = @. uprev + dt*k
  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  if integrator.opts.adaptive
      utilde = @. dt*(integrator.fsalfirst - k)
      integrator.EEst = integrator.opts.internalnorm(
          calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                              integrator.opts.reltol,integrator.opts.internalnorm))
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::MidpointCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::MidpointCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack tmp,k,fsalfirst,atmp = cache
  halfdt = dt/2
  @. tmp = uprev + halfdt*fsalfirst
  f(t+halfdt,tmp,k)
  @. u = uprev + dt*k
  if integrator.opts.adaptive
      @. tmp = dt*(fsalfirst - k)
      calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
                           integrator.opts.reltol,integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  f(t+dt,u,k)
end

function initialize!(integrator,cache::RK4ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::RK4ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  halfdt = dt/2
  k₁ =integrator.fsalfirst
  ttmp = t+halfdt
  k₂ = f(ttmp, @. uprev + halfdt*k₁)
  k₃ = f(ttmp, @. uprev + halfdt*k₂)
  k₄ = f(t+dt, @. uprev + dt*k₃)
  u = @. uprev + (dt/6)*(2*(k₂ + k₃) + (k₁+k₄))
  integrator.fsallast = f(t+dt,u)
  if integrator.opts.adaptive
      # Shampine Solving ODEs and DDEs with Residual Control Estimate
      k₅ = integrator.fsallast
      σ₁ = 1/2 - sqrt(3)/6
      σ₂ = 1/2 + sqrt(3)/6
      p1 = (1-σ₁)*uprev+σ₁*u+σ₁*(σ₁-1)*((1-2σ₁)*(u-uprev)+(σ₁-1)*dt*k₁ + σ₁*dt*k₅)
      p2 = (1-σ₂)*uprev+σ₂*u+σ₂*(σ₂-1)*((1-2σ₂)*(u-uprev)+(σ₂-1)*dt*k₁ + σ₂*dt*k₅)
      pprime1 = k₁ + σ₁*(-4*dt*k₁ - 2*dt*k₅ - 6*uprev +
                σ₁*(3*dt*k₁ + 3*dt*k₅ + 6*uprev - 6*u) + 6*u)/dt
      pprime2 = k₁ + σ₂*(-4*dt*k₁ - 2*dt*k₅ - 6*uprev +
                σ₂*(3*dt*k₁ + 3*dt*k₅ + 6*uprev - 6*u) + 6*u)/dt
      e1 = integrator.opts.internalnorm(calculate_residuals(dt*(f(t+σ₁*dt,p1) - pprime1), uprev, u, integrator.opts.abstol,
      integrator.opts.reltol,integrator.opts.internalnorm))
      e2 = integrator.opts.internalnorm(calculate_residuals(dt*(f(t+σ₂*dt,p2) - pprime2), uprev, u, integrator.opts.abstol, integrator.opts.reltol,
      integrator.opts.internalnorm))
      integrator.EEst = 2.1342*max(e1,e2)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::RK4Cache)
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # pre-start FSAL
end

@muladd function perform_step!(integrator,cache::RK4Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k,atmp = cache
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
  if integrator.opts.adaptive
      # Shampine Solving ODEs and DDEs with Residual Control Estimate
      k₅ = k; p = k₂; pprime = k₃ # Alias some cache arrays
      σ₁ = 1/2 - sqrt(3)/6
      σ₂ = 1/2 + sqrt(3)/6
      @tight_loop_macros for i in eachindex(u)
          @inbounds tmp[i] = (1-σ₁)*uprev[i]+σ₁*u[i]+σ₁*(σ₁-1)*((1-2σ₁)*(u[i]-uprev[i])+(σ₁-1)*dt*k₁[i] + σ₁*dt*k₅[i])
          @inbounds pprime[i] = k₁[i] + σ₁*(-4*dt*k₁[i] - 2*dt*k₅[i] - 6*uprev[i] +
                    σ₁*(3*dt*k₁[i] + 3*dt*k₅[i] + 6*uprev[i] - 6*u[i]) + 6*u[i])/dt
      end
      f(t+σ₁*dt,tmp,p)
      calculate_residuals!(atmp, dt*(p - pprime), uprev, u, integrator.opts.abstol,
                           integrator.opts.reltol,integrator.opts.internalnorm)
      e1 = integrator.opts.internalnorm(atmp)
      @tight_loop_macros for i in eachindex(u)
        @inbounds tmp[i] = (1-σ₂)*uprev[i]+σ₂*u[i]+σ₂*(σ₂-1)*((1-2σ₂)*(u[i]-uprev[i])+(σ₂-1)*dt*k₁[i] + σ₂*dt*k₅[i])
        @inbounds pprime[i] = k₁[i] + σ₂*(-4*dt*k₁[i] - 2*dt*k₅[i] - 6*uprev[i] +
                  σ₂*(3*dt*k₁[i] + 3*dt*k₅[i] + 6*uprev[i] - 6*u[i]) + 6*u[i])/dt
      end
      f(t+σ₂*dt,tmp,p)
      calculate_residuals!(atmp, dt*(p - pprime), uprev, u, integrator.opts.abstol,
                           integrator.opts.reltol,integrator.opts.internalnorm)
      e2 = integrator.opts.internalnorm(atmp)
      integrator.EEst = 2.1342*max(e1,e2)
  end
end


function initialize!(integrator,cache::CarpenterKennedy2N54ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::CarpenterKennedy2N54ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack A2,A3,A4,A5,B1,B2,B3,B4,B5,c2,c3,c4,c5 = cache

  # u1
  tmp = @. dt*integrator.fsalfirst
  u   = @. uprev + B1*tmp
  # u2
  k = f(t+c2*dt, u)
  tmp = @. A2*tmp + dt*k
  u   = @. u + B2*tmp
  # u3
  k = f(t+c3*dt, u)
  tmp = @. A3*tmp + dt*k
  u   = @. u + B3*tmp
  # u4
  k = f(t+c4*dt, u)
  tmp = @. A4*tmp + dt*k
  u   = @. u + B4*tmp
  # u5 = u
  k = f(t+c5*dt, u)
  tmp = @. A5*tmp + dt*k
  u   = @. u + B5*tmp

  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::CarpenterKennedy2N54Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::CarpenterKennedy2N54Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,fsalfirst,tmp,A2,A3,A4,A5,B1,B2,B3,B4,B5,c2,c3,c4,c5 = cache

  # u1
  @. tmp = dt*fsalfirst
  @. u   = uprev + B1*tmp
  # u2
  f(t+c2*dt, u, k)
  @. tmp = A2*tmp + dt*k
  @. u   = u + B2*tmp
  # u3
  f(t+c3*dt, u, k)
  @. tmp = A3*tmp + dt*k
  @. u   = u + B3*tmp
  # u4
  f(t+c4*dt, u, k)
  @. tmp = A4*tmp + dt*k
  @. u   = u + B4*tmp
  # u5 = u
  f(t+c5*dt, u, k)
  @. tmp = A5*tmp + dt*k
  @. u   = u + B5*tmp

  f(t+dt, u, k)
end
