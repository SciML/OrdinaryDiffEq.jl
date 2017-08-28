function initialize!(integrator,cache::SSPRK22ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK22ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator

  # u1 -> stored as u
  u = @. uprev + dt*integrator.fsalfirst
  k = f(t+dt, u)
  # u
  u = @. (uprev + u + dt*k) / 2

  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK22Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK22Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,fsalfirst,stage_limiter!,step_limiter! = cache

  # u1 -> stored as u
  @. u = uprev + dt*fsalfirst
  stage_limiter!(u, f, t+dt)
  f(t+dt, u, k)
  # u
  @. u = (uprev + u + dt*k) / 2
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f(t+dt, u, k)
end


function initialize!(integrator,cache::SSPRK33ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK33ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator

  # u1
  u = @. uprev + dt*integrator.fsalfirst
  k = f(t+dt, u)
  # u2
  u = @. (3*uprev + u + dt*k) / 4
  k = f(t+dt/2, u)
  # u
  u = @. (uprev + 2*u + 2*dt*k) / 3

  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK33Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK33Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,fsalfirst,stage_limiter!,step_limiter! = cache

  # u1
  @. u = uprev + dt*fsalfirst
  stage_limiter!(u, f, t+dt)
  f(t+dt, u, k)
  # u2
  @. u = (3*uprev + u + dt*k) / 4
  stage_limiter!(u, f, t+dt/2)
  f(t+dt/2, u, k)
  # u
  @. u = (uprev + 2*u + 2*dt*k) / 3
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f(t+dt, u, k)
end


function initialize!(integrator,cache::SSPRK53ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK53ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack α30,α32,α40,α43,α52,α54,β10,β21,β32,β43,β54,c1,c2,c3,c4 = cache

  # u1
  tmp = @. uprev + β10 * dt * integrator.fsalfirst
  k = f(t+c1*dt, tmp)
  # u2 -> stored as u
  u = @. tmp + β21 * dt * k
  k = f(t+c2*dt, u)
  # u3
  tmp = @. α30 * uprev + α32 * u + β32 * dt * k
  k = f(t+c3*dt, tmp)
  # u4
  tmp = @. α40 * uprev + α43 * tmp + β43 * dt * k
  k = f(t+c4*dt, tmp)
  # u
  u = @. α52 * u + α54 * tmp + β54 * dt * k

  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK53Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK53Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,tmp,fsalfirst,stage_limiter!,step_limiter!,α30,α32,α40,α43,α52,α54,β10,β21,β32,β43,β54,c1,c2,c3,c4 = cache

  # u1
  @. tmp = uprev + β10 * dt * fsalfirst
  stage_limiter!(tmp, f, t+c1*dt)
  f(t+c1*dt, tmp, k)
  # u2 -> stored as u
  @. u = tmp + β21 * dt * k
  stage_limiter!(u, f, t+c2*dt)
  f(t+c2*dt, u, k)
  # u3
  @. tmp = α30 * uprev + α32 * u + β32 * dt * k
  stage_limiter!(tmp, f, t+c3*dt)
  f(t+c3*dt, tmp, k)
  # u4
  @. tmp = α40 * uprev + α43 * tmp + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f(t+c4*dt, tmp, k)
  # u
  @. u = α52 * u + α54 * tmp + β54 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f(t+dt, u, k)
end


function initialize!(integrator,cache::SSPRK63ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK63ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack α40,α41,α43,α62,α65,β10,β21,β32,β43,β54,β65,c1,c2,c3,c4,c5 = cache

  # u1 -> stored as u
  u = @. uprev + β10 * dt * integrator.fsalfirst
  k = f(t+c1*dt, u)
  # u2
  u₂ = @. u + β21 * dt * k
  k = f(t+c2*dt, u₂)
  # u3
  tmp = @. u₂ + β32 * dt * k
  k = f(t+c3*dt, tmp)
  # u4
  tmp = @. α40 * uprev + α41 * u + α43 * tmp + β43 * dt * k
  k = f(t+c4*dt, tmp)
  # u5
  tmp = @. tmp + β54 * dt * k
  k = f(t+c5*dt, tmp)
  # u
  u = @. α62 * u₂ + α65 * tmp + β65 * dt * k

  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK63Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK63Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,tmp,u₂,fsalfirst,stage_limiter!,step_limiter!,α40,α41,α43,α62,α65,β10,β21,β32,β43,β54,β65,c1,c2,c3,c4,c5 = cache

  # u1 -> stored as u
  @. u = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u, f, t+c1*dt)
  f(t+c1*dt, u, k)
  # u2
  @. u₂ = u + β21 * dt * k
  stage_limiter!(u₂, f, t+c2*dt)
  f(t+c2*dt, u₂, k)
  # u3
  @. tmp = u₂ + β32 * dt * k
  stage_limiter!(tmp, f, t+c3*dt)
  f(t+c3*dt, tmp, k)
  # u4
  @. tmp = α40 * uprev + α41 * u + α43 * tmp + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f(t+c4*dt, tmp, k)
  # u5
  @. tmp = tmp + β54 * dt * k
  stage_limiter!(tmp, f, t+c5*dt)
  f(t+c5*dt, tmp, k)
  # u
  @. u = α62 * u₂ + α65 * tmp + β65 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f(t+dt, u, k)
end


function initialize!(integrator,cache::SSPRK73ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK73ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack α40,α43,α50,α51,α54,α73,α76,β10,β21,β32,β43,β54,β65,β76,c1,c2,c3,c4,c5,c6 = cache

  # u1
  u₁ = @. uprev + β10 * dt * integrator.fsalfirst
  k = f(t+c1*dt, u₁)
  # u2
  tmp = @. u₁ + β21 * dt * k
  k = f(t+c2*dt, tmp)
  # u3 -> stored as u
  u = @. tmp + β32 * dt * k
  k = f(t+c3*dt, u)
  # u4
  tmp = @. α40 * uprev + α43 * u + β43 * dt * k
  k = f(t+c4*dt, tmp)
  # u5
  tmp = @. α50 * uprev + α51 * u₁ + α54 * tmp + β54 * dt * k
  k = f(t+c5*dt, tmp)
  # u6
  tmp = @. tmp + β65 * dt * k
  k = f(t+c6*dt, tmp)
  # u
  u = @. α73 * u + α76 * tmp + β76 * dt * k

  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK73Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK73Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,tmp,u₁,fsalfirst,stage_limiter!,step_limiter!,α40,α43,α50,α51,α54,α73,α76,β10,β21,β32,β43,β54,β65,β76,c1,c2,c3,c4,c5,c6 = cache

  # u1
  @. u₁ = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u₁, f, t+c1*dt)
  f(t+c1*dt, u₁, k)
  # u2
  @. tmp = u₁ + β21 * dt * k
  stage_limiter!(tmp, f, t+c2*dt)
  f(t+c2*dt, tmp, k)
  # u3 -> stored as u
  @. u = tmp + β32 * dt * k
  stage_limiter!(u, f, t+c3*dt)
  f(t+c3*dt, u, k)
  # u4
  @. tmp = α40 * uprev + α43 * u + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f(t+c4*dt, tmp, k)
  # u5
  @. tmp = α50 * uprev + α51 * u₁ + α54 * tmp + β54 * dt * k
  stage_limiter!(tmp, f, t+c5*dt)
  f(t+c5*dt, tmp, k)
  # u6
  @. tmp = tmp + β65 * dt * k
  stage_limiter!(tmp, f, t+c6*dt)
  f(t+c6*dt, tmp, k)
  # u
  @. u = α73 * u + α76 * tmp + β76 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f(t+dt, u, k)
end


function initialize!(integrator,cache::SSPRK83ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK83ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack α50,α51,α54,α61,α65,α72,α73,α76,β10,β21,β32,β43,β54,β65,β76,β87,c1,c2,c3,c4,c5,c6,c7 = cache

  # u1 -> save as u
  u = @. uprev + β10 * dt * integrator.fsalfirst
  k = f(t+c1*dt, u)
  # u2
  u₂ = @. u + β21 * dt * k
  k = f(t+c2*dt, u₂)
  # u3
  u₃ = @. u₂ + β32 * dt * k
  k = f(t+c3*dt, u₃)
  # u4
  tmp = @. u₃ + β43 * dt * k
  k = f(t+c4*dt, tmp)
  # u5
  tmp = @. α50 * uprev + α51 * u + α54 * tmp + β54 * dt * k
  k = f(t+c5*dt, tmp)
  # u6
  tmp = @. α61 * u + α65 * tmp + β65 * dt * k
  k = f(t+c6*dt, tmp)
  # u7
  tmp = @. α72 * u₂ + α73 * u₃ + α76 * tmp + β76 * dt * k
  k = f(t+c7*dt, tmp)
  # u
  u = @. tmp + β87 * dt * k

  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK83Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK83Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,tmp,u₂,u₃,fsalfirst,stage_limiter!,step_limiter!,α50,α51,α54,α61,α65,α72,α73,α76,β10,β21,β32,β43,β54,β65,β76,β87,c1,c2,c3,c4,c5,c6,c7 = cache

  # u1 -> save as u
  @. u = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u, f, t+c1*dt)
  f(t+c1*dt, u, k)
  # u2
  @. u₂ = u + β21 * dt * k
  stage_limiter!(u₂, f, t+c2*dt)
  f(t+c2*dt, u₂, k)
  # u3
  @. u₃ = u₂ + β32 * dt * k
  stage_limiter!(u₃, f, t+c3*dt)
  f(t+c3*dt, u₃, k)
  # u4
  @. tmp = u₃ + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f(t+c4*dt, tmp, k)
  # u5
  @. tmp = α50 * uprev + α51 * u + α54 * tmp + β54 * dt * k
  stage_limiter!(tmp, f, t+c5*dt)
  f(t+c5*dt, tmp, k)
  # u6
  @. tmp = α61 * u + α65 * tmp + β65 * dt * k
  stage_limiter!(tmp, f, t+c6*dt)
  f(t+c6*dt, tmp, k)
  # u7
  @. tmp = α72 * u₂ + α73 * u₃ + α76 * tmp + β76 * dt * k
  stage_limiter!(tmp, f, t+c7*dt)
  f(t+c7*dt, tmp, k)
  # u
  @. u = tmp + β87 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f(t+dt, u, k)
end


function initialize!(integrator,cache::SSPRK432ConstantCache)
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK432ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  dt_2 = dt / 2
  tmp = @. uprev + dt_2*integrator.fsalfirst # u1
  k = f(t+dt_2, tmp)
  tmp = @. tmp   + dt_2*k # u2
  k = f(t+dt, tmp)
  tmp = @. tmp   + dt_2*k
  if integrator.opts.adaptive
    utilde = @. (uprev + 2*tmp) / 3
  end
  tmp = @. (2*uprev + tmp) / 3 #u3
  k = f(t+dt_2, tmp)
  u = @. tmp + dt_2*k
  integrator.fsallast = f(t+dt,u)
  if integrator.opts.adaptive
    integrator.EEst = integrator.opts.internalnorm(@. (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK432Cache)
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # Pre-start fsal
end

@muladd function perform_step!(integrator,cache::SSPRK432Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,tmp,fsalfirst,utilde,atmp,stage_limiter!,step_limiter! = cache
  dt_2 = dt / 2

  # u1
  @. tmp = uprev + dt_2*fsalfirst
  stage_limiter!(tmp, f, t+dt_2)
  f(t+dt_2, tmp, k)
  # u2
  @. tmp = tmp + dt_2*k
  stage_limiter!(tmp, f, t+dt)
  f(t+dt, tmp, k)
  #
  @. tmp = tmp + dt_2*k
  stage_limiter!(tmp, f, t+dt+dt_2)
  if integrator.opts.adaptive
    @. utilde = (uprev + 2*tmp) / 3
  end
  # u3
  @. tmp = (2*uprev + tmp) / 3
  f(t+dt_2, tmp, k)
  #
  @. u = tmp + dt_2*k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)

  if integrator.opts.adaptive
    @. atmp = (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  f(t+dt, u, k)
end


function initialize!(integrator,cache::SSPRK932ConstantCache)
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK932ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  dt_6 = dt / 6
  dt_3 = dt / 3
  dt_2 = dt / 2

  tmp = @. uprev + dt_6*integrator.fsalfirst # u1
  k = f(t+dt_6, tmp)
  tmp = @. tmp   + dt_6*k # u2
  k = f(t+dt_3, tmp)
  tmp = @. tmp   + dt_6*k # u3
  k = f(t+dt_2, tmp)
  tmp = @. tmp   + dt_6*k # u4
  k = f(t+2*dt_3, tmp)
  tmp = @. tmp   + dt_6*k # u5
  k = f(t+5*dt_6, tmp)
  tmp = @. tmp   + dt_6*k # u6
  if integrator.opts.adaptive
    k = f(t+dt, tmp)
    utilde = @. (uprev + 6*tmp + 6*dt*k) / 7
  end
  tmp = @. (3*uprev + dt_2*integrator.fsalfirst + 2*tmp) / 5 # u6*
  k = f(t+dt_2, tmp)
  tmp = @. tmp   + dt_6*k # u7*
  k = f(t+2*dt_3, tmp)
  tmp = @. tmp   + dt_6*k # u8*
  k = f(t+5*dt_6, tmp)
  u = @. tmp     + dt_6*k # u9*
  integrator.fsallast = f(t+dt,u)
  if integrator.opts.adaptive
    integrator.EEst = integrator.opts.internalnorm(@. (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK932Cache)
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # Pre-start fsal
end

@muladd function perform_step!(integrator,cache::SSPRK932Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,tmp,fsalfirst,utilde,atmp,stage_limiter!,step_limiter! = cache
  dt_6 = dt / 6
  dt_3 = dt / 3
  dt_2 = dt / 2

  # u1
  @. tmp = uprev + dt_6*fsalfirst
  stage_limiter!(tmp, f, t+dt_6)
  f(t+dt_6, tmp, k)
  # u2
  @. tmp = tmp + dt_6*k
  stage_limiter!(tmp, f, t+dt_3)
  f(t+dt_3, tmp, k)
  # u3
  @. tmp = tmp + dt_6*k
  stage_limiter!(tmp, f, t+dt_2)
  f(t+dt_2, tmp, k)
  # u4
  @. tmp = tmp + dt_6*k
  stage_limiter!(tmp, f, t+2*dt_3)
  f(t+2*dt_3, tmp, k)
  # u5
  @. tmp = tmp + dt_6*k
  stage_limiter!(tmp, f, t+5*dt_6)
  f(t+5*dt_6, tmp, k)
  # u6
  @. tmp = tmp + dt_6*k
  if integrator.opts.adaptive
    stage_limiter!(tmp, f, t+dt)
    f(t+dt, tmp, k)
    @. utilde = (uprev + 6*tmp + 6*dt*k) / 7
    stage_limiter!(utilde, f, t+dt)
  end
  # u6*
  @. tmp = (3*uprev + dt_2*fsalfirst + 2*tmp) / 5
  stage_limiter!(tmp, f, t+dt_6)
  f(t+dt_2, tmp, k)
  # u7*
  @. tmp = tmp + dt_6*k
  stage_limiter!(tmp, f, t+2*dt_3)
  f(t+2*dt_3, tmp, k)
  # u8*
  @. tmp = tmp + dt_6*k
  stage_limiter!(tmp, f, t+5*dt_6)
  f(t+5*dt_6, tmp, k)
  # u9*
  @. u = tmp + dt_6*k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)

  if integrator.opts.adaptive
    @. atmp = (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  f(t+dt, u, k)
end


function initialize!(integrator,cache::SSPRK54ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::SSPRK54ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack β10,α20,α21,β21,α30,α32,β32,α40,α43,β43,α52,α53,β53,α54,β54,c1,c2,c3,c4 = cache

  # u₁
  u₂ = @. uprev + β10 * dt * integrator.fsalfirst
  k = f(t+c1*dt, u₂)
  # u₂
  u₂ = @. α20 * uprev + α21 * u₂ + β21 * dt * k
  k = f(t+c2*dt, u₂)
  # u₃
  u₃ = @. α30 * uprev + α32 * u₂ + β32 * dt * k
  k₃ = f(t+c3*dt, u₃)
  # u₄
  u₄ = @. α40 * uprev + α43 * u₃ + β43 * dt * k₃
  k = f(t+c4*dt, u₄)
  # u
  u = @. α52 * u₂ + α53 * u₃ + β53 * dt * k₃ + α54 * u₄ + β54 * dt * k

  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK54Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK54Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,k₃,u₂,u₃,u₄,fsalfirst,stage_limiter!,step_limiter!,β10,α20,α21,β21,α30,α32,β32,α40,α43,β43,α52,α53,β53,α54,β54,c1,c2,c3,c4 = cache

  # u₁
  @. u₂ = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u₂, f, t+c1*dt)
  f(t+c1*dt, u₂, k)
  # u₂
  @. u₂ = α20 * uprev + α21 * u₂ + β21 * dt * k
  stage_limiter!(u₂, f, t+c2*dt)
  f(t+c2*dt, u₂, k)
  # u₃
  @. u₃ = α30 * uprev + α32 * u₂ + β32 * dt * k
  stage_limiter!(u₃, f, t+c3*dt)
  f(t+c3*dt, u₃, k₃)
  # u₄
  @. u₄ = α40 * uprev + α43 * u₃ + β43 * dt * k₃
  stage_limiter!(u₄, f, t+c4*dt)
  f(t+c4*dt, u₄, k)
  # u
  @. u = α52 * u₂ + α53 * u₃ + β53 * dt * k₃ + α54 * u₄ + β54 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f(t+dt, u, k)
end


function initialize!(integrator,cache::SSPRK104ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::SSPRK104ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  dt_6 = dt/6
  dt_3 = dt/3
  dt_2 = dt/2

  tmp = @. uprev + dt_6 * integrator.fsalfirst # u₁
  k = f(t+dt_6, tmp)
  tmp = @. tmp + dt_6 * k # u₂
  k = f(t+dt_3, tmp)
  tmp = @. tmp + dt_6 * k # u₃
  k = f(t+dt_2, tmp)
  u = @. tmp + dt_6 * k # u₄
  k₄ = f(t+2*dt_3, u)
  tmp = @. (3*uprev + 2*u + dt_3 * k₄) / 5 # u₅
  k = f(t+dt_3, tmp)
  tmp = @. tmp + dt_6 * k # u₆
  k = f(t+dt_2, tmp)
  tmp = @. tmp + dt_6 * k # u₇
  k = f(t+2*dt_3, tmp)
  tmp = @. tmp + dt_6 * k # u₈
  k = f(t+5*dt_6, tmp)
  tmp = @. tmp + dt_6 * k # u₉
  k = f(t+dt, tmp)
  u = @. (uprev + 9*(u + dt_6*k₄) + 15*(tmp + dt_6*k)) / 25

  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK104Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK104Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,k₄,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  dt_6 = dt/6
  dt_3 = dt/3
  dt_2 = dt/2

  @. tmp = uprev + dt_6 * integrator.fsalfirst
  stage_limiter!(tmp, f, t+dt_6)
  f(t+dt_6, tmp, k)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt_3)
  f(t+dt_3, tmp, k)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt_2)
  f(t+dt_2, tmp, k)
  @. u = tmp + dt_6 * k
  stage_limiter!(u, f, t+2*dt_3)
  f(t+2*dt_3, u, k₄)
  @. tmp = (3*uprev + 2*u + 2*dt_6 * k₄) / 5
  stage_limiter!(tmp, f, t+dt_3)
  f(t+dt_3, tmp, k)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt_2)
  f(t+dt_2, tmp, k)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+2*dt_3)
  f(t+2*dt_3, tmp, k)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+5*dt_6)
  f(t+5*dt_6, tmp, k)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt)
  f(t+dt, tmp, k)
  @. u = (uprev + 9*(u + dt_6*k₄) + 15*(tmp + dt_6*k)) / 25
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f(t+dt,u,k)
end
