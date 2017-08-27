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
  tmp = @. uprev + dt*integrator.fsalfirst
  k = f(t+dt,tmp)
  u = @. (uprev + tmp + dt*k) / 2
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
  @unpack k,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  @. tmp = uprev + dt*integrator.fsalfirst
  stage_limiter!(tmp, f, t+dt)
  f(t+dt,tmp,k)
  @. u = (uprev + tmp + dt*k) / 2
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f(t+dt,u,k)
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
  tmp = @. uprev + dt*integrator.fsalfirst
  k = f(t+dt,tmp)
  tmp = @. (3*uprev + tmp + dt*k) / 4
  k = f(t+dt/2,tmp)
  u = @. (uprev + 2*tmp + 2*dt*k) / 3
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
  @unpack k,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  @. tmp = uprev + dt*fsalfirst
  stage_limiter!(tmp, f, t+dt)
  f(t+dt,tmp,k)
  @. tmp = (3*uprev + tmp + dt*k) / 4
  stage_limiter!(tmp, f, t+dt/2)
  f(t+dt/2,tmp,k)
  @. u = (uprev + 2*tmp + 2*dt*k) / 3
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  f(t+dt,u,k)
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
