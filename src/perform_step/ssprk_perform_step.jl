function initialize!(integrator,cache::SSPRK22ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK22ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator

  # u1 -> stored as u
  u = uprev + dt*integrator.fsalfirst
  k = f(u, p, t+dt)
  # u
  u = (uprev + u + dt*k) / 2

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK22Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK22Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,stage_limiter!,step_limiter! = cache

  # u1 -> stored as u
  @. u = uprev + dt*fsalfirst
  stage_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
  # u
  @. u = (uprev + u + dt*k) / 2
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK33ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK33ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator

  # u1
  u = uprev + dt*integrator.fsalfirst
  k = f(u, p, t+dt)
  # u2
  u = (3*uprev + u + dt*k) / 4
  k = f(u,p,t+dt/2)
  # u
  u = (uprev + 2*u + 2*dt*k) / 3

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK33Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK33Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,stage_limiter!,step_limiter! = cache

  # u1
  @. u = uprev + dt*fsalfirst
  stage_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
  # u2
  @. u = (3*uprev + u + dt*k) / 4
  stage_limiter!(u, f, t+dt/2)
  f(k,u,p,t+dt/2)
  # u
  @. u = (uprev + 2*u + 2*dt*k) / 3
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK53ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK53ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α30,α32,α40,α43,α52,α54,β10,β21,β32,β43,β54,c1,c2,c3,c4 = cache

  # u1
  tmp = uprev + β10 * dt * integrator.fsalfirst
  k = f(tmp, p, t+c1*dt)
  # u2 -> stored as u
  u = tmp + β21 * dt * k
  k = f(u, p, t+c2*dt)
  # u3
  tmp = α30 * uprev + α32 * u + β32 * dt * k
  k = f(tmp, p, t+c3*dt)
  # u4
  tmp = α40 * uprev + α43 * tmp + β43 * dt * k
  k = f(tmp, p, t+c4*dt)
  # u
  u = α52 * u + α54 * tmp + β54 * dt * k

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK53Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK53Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack α30,α32,α40,α43,α52,α54,β10,β21,β32,β43,β54,c1,c2,c3,c4 = cache.tab

  # u1
  @. tmp = uprev + β10 * dt * fsalfirst
  stage_limiter!(tmp, f, t+c1*dt)
  f( k,  tmp, p, t+c1*dt)
  # u2 -> stored as u
  @. u = tmp + β21 * dt * k
  stage_limiter!(u, f, t+c2*dt)
  f( k,  u, p, t+c2*dt)
  # u3
  @. tmp = α30 * uprev + α32 * u + β32 * dt * k
  stage_limiter!(tmp, f, t+c3*dt)
  f( k,  tmp, p, t+c3*dt)
  # u4
  @. tmp = α40 * uprev + α43 * tmp + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f( k,  tmp, p, t+c4*dt)
  # u
  @. u = α52 * u + α54 * tmp + β54 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK63ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK63ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α40,α41,α43,α62,α65,β10,β21,β32,β43,β54,β65,c1,c2,c3,c4,c5 = cache

  # u1 -> stored as u
  u = uprev + β10 * dt * integrator.fsalfirst
  k = f(u, p, t+c1*dt)
  # u2
  u₂ = u + β21 * dt * k
  k = f(u₂,p,t+c2*dt)
  # u3
  tmp = u₂ + β32 * dt * k
  k = f(tmp, p, t+c3*dt)
  # u4
  tmp = α40 * uprev + α41 * u + α43 * tmp + β43 * dt * k
  k = f(tmp, p, t+c4*dt)
  # u5
  tmp = tmp + β54 * dt * k
  k = f(tmp, p, t+c5*dt)
  # u
  u = α62 * u₂ + α65 * tmp + β65 * dt * k

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK63Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK63Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,u₂,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack α40,α41,α43,α62,α65,β10,β21,β32,β43,β54,β65,c1,c2,c3,c4,c5 = cache.tab

  # u1 -> stored as u
  @. u = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u, f, t+c1*dt)
  f( k,  u, p, t+c1*dt)
  # u2
  @. u₂ = u + β21 * dt * k
  stage_limiter!(u₂, f, t+c2*dt)
  f(k,u₂,p,t+c2*dt)
  # u3
  @. tmp = u₂ + β32 * dt * k
  stage_limiter!(tmp, f, t+c3*dt)
  f( k,  tmp, p, t+c3*dt)
  # u4
  @. tmp = α40 * uprev + α41 * u + α43 * tmp + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f( k,  tmp, p, t+c4*dt)
  # u5
  @. tmp = tmp + β54 * dt * k
  stage_limiter!(tmp, f, t+c5*dt)
  f( k,  tmp, p, t+c5*dt)
  # u
  @. u = α62 * u₂ + α65 * tmp + β65 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK73ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK73ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α40,α43,α50,α51,α54,α73,α76,β10,β21,β32,β43,β54,β65,β76,c1,c2,c3,c4,c5,c6 = cache

  # u1
  u₁ = uprev + β10 * dt * integrator.fsalfirst
  k = f(u₁,p,t+c1*dt)
  # u2
  tmp = u₁ + β21 * dt * k
  k = f(tmp, p, t+c2*dt)
  # u3 -> stored as u
  u = tmp + β32 * dt * k
  k = f(u, p, t+c3*dt)
  # u4
  tmp = α40 * uprev + α43 * u + β43 * dt * k
  k = f(tmp, p, t+c4*dt)
  # u5
  tmp = α50 * uprev + α51 * u₁ + α54 * tmp + β54 * dt * k
  k = f(tmp, p, t+c5*dt)
  # u6
  tmp = tmp + β65 * dt * k
  k = f(tmp, p, t+c6*dt)
  # u
  u = α73 * u + α76 * tmp + β76 * dt * k

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK73Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK73Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,u₁,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack α40,α43,α50,α51,α54,α73,α76,β10,β21,β32,β43,β54,β65,β76,c1,c2,c3,c4,c5,c6 = cache.tab

  # u1
  @. u₁ = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u₁, f, t+c1*dt)
  f(k,u₁,p,t+c1*dt)
  # u2
  @. tmp = u₁ + β21 * dt * k
  stage_limiter!(tmp, f, t+c2*dt)
  f( k,  tmp, p, t+c2*dt)
  # u3 -> stored as u
  @. u = tmp + β32 * dt * k
  stage_limiter!(u, f, t+c3*dt)
  f( k,  u, p, t+c3*dt)
  # u4
  @. tmp = α40 * uprev + α43 * u + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f( k,  tmp, p, t+c4*dt)
  # u5
  @. tmp = α50 * uprev + α51 * u₁ + α54 * tmp + β54 * dt * k
  stage_limiter!(tmp, f, t+c5*dt)
  f( k,  tmp, p, t+c5*dt)
  # u6
  @. tmp = tmp + β65 * dt * k
  stage_limiter!(tmp, f, t+c6*dt)
  f( k,  tmp, p, t+c6*dt)
  # u
  @. u = α73 * u + α76 * tmp + β76 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK83ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK83ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α50,α51,α54,α61,α65,α72,α73,α76,β10,β21,β32,β43,β54,β65,β76,β87,c1,c2,c3,c4,c5,c6,c7 = cache

  # u1 -> save as u
  u = uprev + β10 * dt * integrator.fsalfirst
  k = f(u, p, t+c1*dt)
  # u2
  u₂ = u + β21 * dt * k
  k = f(u₂,p,t+c2*dt)
  # u3
  u₃ = u₂ + β32 * dt * k
  k = f(u₃,p,t+c3*dt)
  # u4
  tmp = u₃ + β43 * dt * k
  k = f(tmp, p, t+c4*dt)
  # u5
  tmp = α50 * uprev + α51 * u + α54 * tmp + β54 * dt * k
  k = f(tmp, p, t+c5*dt)
  # u6
  tmp = α61 * u + α65 * tmp + β65 * dt * k
  k = f(tmp, p, t+c6*dt)
  # u7
  tmp = α72 * u₂ + α73 * u₃ + α76 * tmp + β76 * dt * k
  k = f(tmp, p, t+c7*dt)
  # u
  u = tmp + β87 * dt * k

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK83Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK83Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,u₂,u₃,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack α50,α51,α54,α61,α65,α72,α73,α76,β10,β21,β32,β43,β54,β65,β76,β87,c1,c2,c3,c4,c5,c6,c7 = cache.tab

  # u1 -> save as u
  @. u = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u, f, t+c1*dt)
  f( k,  u, p, t+c1*dt)
  # u2
  @. u₂ = u + β21 * dt * k
  stage_limiter!(u₂, f, t+c2*dt)
  f(k,u₂,p,t+c2*dt)
  # u3
  @. u₃ = u₂ + β32 * dt * k
  stage_limiter!(u₃, f, t+c3*dt)
  f(k,u₃,p,t+c3*dt)
  # u4
  @. tmp = u₃ + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f( k,  tmp, p, t+c4*dt)
  # u5
  @. tmp = α50 * uprev + α51 * u + α54 * tmp + β54 * dt * k
  stage_limiter!(tmp, f, t+c5*dt)
  f( k,  tmp, p, t+c5*dt)
  # u6
  @. tmp = α61 * u + α65 * tmp + β65 * dt * k
  stage_limiter!(tmp, f, t+c6*dt)
  f( k,  tmp, p, t+c6*dt)
  # u7
  @. tmp = α72 * u₂ + α73 * u₃ + α76 * tmp + β76 * dt * k
  stage_limiter!(tmp, f, t+c7*dt)
  f( k,  tmp, p, t+c7*dt)
  # u
  @. u = tmp + β87 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK432ConstantCache)
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK432ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  dt_2 = dt / 2

  # u1
  u = uprev + dt_2*integrator.fsalfirst
  k = f(u, p, t+dt_2)
  # u2
  u = u + dt_2*k
  k = f(u, p, t+dt)
  u = u + dt_2*k
  if integrator.opts.adaptive
    utilde = (uprev + 2*u) / 3
  end
  # u3
  u = (2*uprev + u) / 3
  k = f(u, p, t+dt_2)
  # u
  u = u + dt_2*k

  integrator.fsallast = f(u, p, t+dt)
  if integrator.opts.adaptive
    utilde = utilde - u
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK432Cache)
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
end

@muladd function perform_step!(integrator,cache::SSPRK432Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,utilde,atmp,stage_limiter!,step_limiter! = cache
  dt_2 = dt / 2

  # u1
  @. u = uprev + dt_2*fsalfirst
  stage_limiter!(u, f, t+dt_2)
  f( k,  u, p, t+dt_2)
  # u2
  @. u = u + dt_2*k
  stage_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
  #
  @. u = u + dt_2*k
  stage_limiter!(u, f, t+dt+dt_2)
  if integrator.opts.adaptive
    @. utilde = (uprev + 2*u) / 3
  end
  # u3
  @. u = (2*uprev + u) / 3
  f( k,  u, p, t+dt_2)
  #
  @. u = u + dt_2*k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)

  if integrator.opts.adaptive
    @. utilde = utilde - u
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK932ConstantCache)
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK932ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  dt_6 = dt / 6
  dt_3 = dt / 3
  dt_2 = dt / 2

  # u1
  u = uprev + dt_6*integrator.fsalfirst
  k = f(u, p, t+dt_6)
  # u2
  u = u + dt_6*k
  k = f(u, p, t+dt_3)
  # u3
  u = u + dt_6*k
  k = f(u, p, t+dt_2)
  # u4
  u = u + dt_6*k
  k = f(u, p, t+2*dt_3)
  # u5
  u = u + dt_6*k
  k = f(u, p, t+5*dt_6)
  # u6
  u = u + dt_6*k
  if integrator.opts.adaptive
    k = f(u, p, t+dt)
    utilde = (uprev + 6*u + 6*dt*k) / 7
  end
  # u6*
  u = (3*uprev + dt_2*integrator.fsalfirst + 2*u) / 5
  k = f(u, p, t+dt_2)
  # u7*
  u = u + dt_6*k
  k = f(u, p, t+2*dt_3)
  # u8*
  u = u + dt_6*k
  k = f(u, p, t+5*dt_6)
  # u
  u = u + dt_6*k

  integrator.fsallast = f(u, p, t+dt)
  if integrator.opts.adaptive
    utilde = utilde - u
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK932Cache)
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
end

@muladd function perform_step!(integrator,cache::SSPRK932Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,utilde,atmp,stage_limiter!,step_limiter! = cache
  dt_6 = dt / 6
  dt_3 = dt / 3
  dt_2 = dt / 2

  # u1
  @. u = uprev + dt_6*fsalfirst
  stage_limiter!(u, f, t+dt_6)
  f( k,  u, p, t+dt_6)
  # u2
  @. u = u + dt_6*k
  stage_limiter!(u, f, t+dt_3)
  f( k,  u, p, t+dt_3)
  # u3
  @. u = u + dt_6*k
  stage_limiter!(u, f, t+dt_2)
  f( k,  u, p, t+dt_2)
  # u4
  @. u = u + dt_6*k
  stage_limiter!(u, f, t+2*dt_3)
  f( k,  u, p, t+2*dt_3)
  # u5
  @. u = u + dt_6*k
  stage_limiter!(u, f, t+5*dt_6)
  f( k,  u, p, t+5*dt_6)
  # u6
  @. u = u + dt_6*k
  if integrator.opts.adaptive
    stage_limiter!(u, f, t+dt)
    f( k,  u, p, t+dt)
    @. utilde = (uprev + 6*u + 6*dt*k) / 7
    stage_limiter!(utilde, f, t+dt)
  end
  # u6*
  @. u = (3*uprev + dt_2*fsalfirst + 2*u) / 5
  stage_limiter!(u, f, t+dt_6)
  f( k,  u, p, t+dt_2)
  # u7*
  @. u = u + dt_6*k
  stage_limiter!(u, f, t+2*dt_3)
  f( k,  u, p, t+2*dt_3)
  # u8*
  @. u = u + dt_6*k
  stage_limiter!(u, f, t+5*dt_6)
  f( k,  u, p, t+5*dt_6)
  # u9*
  @. u = u + dt_6*k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)

  if integrator.opts.adaptive
    @. utilde = utilde - u
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK54ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::SSPRK54ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack β10,α20,α21,β21,α30,α32,β32,α40,α43,β43,α52,α53,β53,α54,β54,c1,c2,c3,c4 = cache

  # u₁
  u₂ = uprev + β10 * dt * integrator.fsalfirst
  k = f(u₂,p,t+c1*dt)
  # u₂
  u₂ = α20 * uprev + α21 * u₂ + β21 * dt * k
  k = f(u₂,p,t+c2*dt)
  # u₃
  u₃ = α30 * uprev + α32 * u₂ + β32 * dt * k
  k₃ = f(u₃,p,t+c3*dt)
  # u₄ -> stored as tmp
  tmp = α40 * uprev + α43 * u₃ + β43 * dt * k₃
  k = f(tmp, p, t+c4*dt)
  # u
  u = α52 * u₂ + α53 * u₃ + β53 * dt * k₃ + α54 * tmp + β54 * dt * k

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK54Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK54Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,k₃,u₂,u₃,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack β10,α20,α21,β21,α30,α32,β32,α40,α43,β43,α52,α53,β53,α54,β54,c1,c2,c3,c4 = cache.tab

  # u₁
  @. u₂ = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u₂, f, t+c1*dt)
  f(k,u₂,p,t+c1*dt)
  # u₂
  @. u₂ = α20 * uprev + α21 * u₂ + β21 * dt * k
  stage_limiter!(u₂, f, t+c2*dt)
  f(k,u₂,p,t+c2*dt)
  # u₃
  @. u₃ = α30 * uprev + α32 * u₂ + β32 * dt * k
  stage_limiter!(u₃, f, t+c3*dt)
  f(k₃,u₃,p,t+c3*dt)
  # u₄ -> stored as tmp
  @. tmp = α40 * uprev + α43 * u₃ + β43 * dt * k₃
  stage_limiter!(tmp, f, t+c4*dt)
  f( k,  tmp, p, t+c4*dt)
  # u
  @. u = α52 * u₂ + α53 * u₃ + β53 * dt * k₃ + α54 * tmp + β54 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK104ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::SSPRK104ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  dt_6 = dt/6
  dt_3 = dt/3
  dt_2 = dt/2

  tmp = uprev + dt_6 * integrator.fsalfirst # u₁
  k = f(tmp, p, t+dt_6)
  tmp = tmp + dt_6 * k # u₂
  k = f(tmp, p, t+dt_3)
  tmp = tmp + dt_6 * k # u₃
  k = f(tmp, p, t+dt_2)
  u = tmp + dt_6 * k # u₄
  k₄ = f(u, p, t+2*dt_3)
  tmp = (3*uprev + 2*u + dt_3 * k₄) / 5 # u₅
  k = f(tmp, p, t+dt_3)
  tmp = tmp + dt_6 * k # u₆
  k = f(tmp, p, t+dt_2)
  tmp = tmp + dt_6 * k # u₇
  k = f(tmp, p, t+2*dt_3)
  tmp = tmp + dt_6 * k # u₈
  k = f(tmp, p, t+5*dt_6)
  tmp = tmp + dt_6 * k # u₉
  k = f(tmp, p, t+dt)
  u = (uprev + 9*(u + dt_6*k₄) + 15*(tmp + dt_6*k)) / 25

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK104Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK104Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,k₄,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  dt_6 = dt/6
  dt_3 = dt/3
  dt_2 = dt/2

  @. tmp = uprev + dt_6 * integrator.fsalfirst
  stage_limiter!(tmp, f, t+dt_6)
  f( k,  tmp, p, t+dt_6)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt_3)
  f( k,  tmp, p, t+dt_3)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt_2)
  f( k,  tmp, p, t+dt_2)
  @. u = tmp + dt_6 * k
  stage_limiter!(u, f, t+2*dt_3)
  f(k₄,u,p,t+2*dt_3)
  @. tmp = (3*uprev + 2*u + 2*dt_6 * k₄) / 5
  stage_limiter!(tmp, f, t+dt_3)
  f( k,  tmp, p, t+dt_3)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt_2)
  f( k,  tmp, p, t+dt_2)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+2*dt_3)
  f( k,  tmp, p, t+2*dt_3)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+5*dt_6)
  f( k,  tmp, p, t+5*dt_6)
  @. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt)
  f( k,  tmp, p, t+dt)
  @. u = (uprev + 9*(u + dt_6*k₄) + 15*(tmp + dt_6*k)) / 25
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f(k, u, p, t+dt)
end
