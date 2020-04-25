function initialize!(integrator,cache::SSPRK22ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
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
  integrator.destats.nf += 2
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK22Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,stage_limiter!,step_limiter! = cache

  # u1 -> stored as u
  @.. u = uprev + dt*fsalfirst
  stage_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
  # u
  @.. u = (uprev + u + dt*k) / 2
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 2
  f( k,  u, p, t+dt)
end

function initialize!(integrator,cache::KYKSSPRK42Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::KYKSSPRK42Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp, fsalfirst = cache
  @unpack α20, α21, α30, α32, α40, α43, β10, β21, β30, β32, β40, β43, c1, c2, c3 = cache.tab

  δ = fsalfirst
  # u1 -> stored as u
  @.. u = uprev + dt*β10*δ
  f(k, u, p, t+c1*dt)
  # u2
  @.. tmp = α20*uprev + α21*u + dt*β21*k
  f(k, tmp, p, t+c2*dt)
  # u3
  @.. tmp = α30*uprev + α32*tmp + dt*β30*δ + dt*β32*k
  f(k, tmp, p, t+c3*dt)
  # u
  @.. u = α40*uprev + α43*tmp + dt*β40*δ + dt*β43*k
  f(k, u, p, t+dt)
end

function initialize!(integrator,cache::KYKSSPRK42ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::KYKSSPRK42ConstantCache)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α20, α21, α30, α32, α40, α43, β10, β21, β30, β32, β40, β43, c1, c2, c3 = cache

  #u1
  δ = integrator.fsalfirst
  u = uprev + dt*β10*δ
  k = f(u, p, t + c1*dt)
  #u2
  tmp = α20*uprev + α21*u + dt*β21*k
  k = f(tmp, p, t + c2*dt)
  #u3
  tmp = α30*uprev + α32*tmp + dt*β30*δ + dt*β32*k
  k = f(tmp, p, t + c3*dt)
  #u
  u = α40*uprev + α43*tmp + dt*β40*δ + dt*β43*k

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SHLDDRK52ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::SHLDDRK52ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α2,α3,α4,α5,β1,β2,β3,β4,β5,c2,c3,c4,c5 = cache

  # u1
  tmp = dt*integrator.fsalfirst
  u   = uprev + β1*tmp
  # u2
  tmp = α2*tmp + dt*f(u, p, t+c2*dt)
  u   = u + β2*tmp
  # u3
  tmp = α3*tmp + dt*f(u, p, t+c3*dt)
  u   = u + β3*tmp
  # u4
  tmp = α4*tmp + dt*f(u, p, t+c4*dt)
  u   = u + β4*tmp
  # u5 = u
  tmp = α5*tmp + dt*f(u, p, t+c5*dt)
  u   = u + β5*tmp

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::SHLDDRK52Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SHLDDRK52Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,tmp = cache
  @unpack α2,α3,α4,α5,β1,β2,β3,β4,β5,c2,c3,c4,c5 = cache.tab

  # u1
  @. tmp = dt*fsalfirst
  @. u   = uprev + β1*tmp
  # u2
  f( k,  u, p, t+c2*dt)
  @. tmp = α2*tmp + dt*k
  @. u   = u + β2*tmp
  # u3
  f( k,  u, p, t+c3*dt)
  @. tmp = α3*tmp + dt*k
  @. u   = u + β3*tmp
  # u4
  f( k,  u, p, t+c4*dt)
  @. tmp = α4*tmp + dt*k
  @. u   = u + β4*tmp
  # u5 = u
  f( k,  u, p, t+c5*dt)
  @. tmp = α5*tmp + dt*k
  @. u   = u + β5*tmp

  f( k,  u, p, t+dt)
end

function initialize!(integrator,cache::SHLDDRK_2NConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::SHLDDRK_2NConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α21,α31,α41,α51,β11,β21,β31,β41,β51,c21,c31,c41,c51,α22,α32,α42,α52,α62,β12,β22,β32,β42,β52,β62,c22,c32,c42,c52,c62= cache

  if integrator.u_modified
    cache.step = 1
  end
  # cnt = cache.step

  if cache.step % 2 == 1
    cache.step += 1
    # u1
    tmp = dt*integrator.fsalfirst
    u   = uprev + β11*tmp
    # u2
    tmp = α21*tmp + dt*f(u, p, t+c21*dt)
    u   = u + β21*tmp
    # u3
    tmp = α31*tmp + dt*f(u, p, t+c31*dt)
    u   = u + β31*tmp
    # u4
    tmp = α41*tmp + dt*f(u, p, t+c41*dt)
    u   = u + β41*tmp
    # u5 = u
    tmp = α51*tmp + dt*f(u, p, t+c51*dt)
    u   = u + β51*tmp

  else
    cache.step += 1
    # u1
    tmp = dt*integrator.fsalfirst
    u   = uprev + β12*tmp
    # u2
    tmp = α22*tmp + dt*f(u, p, t+c22*dt)
    u   = u + β22*tmp
    # u3
    tmp = α32*tmp + dt*f(u, p, t+c32*dt)
    u   = u + β32*tmp
    # u4
    tmp = α42*tmp + dt*f(u, p, t+c42*dt)
    u   = u + β42*tmp
    # u5 = u
    tmp = α52*tmp + dt*f(u, p, t+c52*dt)
    u   = u + β52*tmp
    tmp = α62*tmp + dt*f(u, p, t+c62*dt)
    u   = u + β62*tmp
  end

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::SHLDDRK_2NCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SHLDDRK_2NCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,tmp = cache
  @unpack α21,α31,α41,α51,β11,β21,β31,β41,β51,c21,c31,c41,c51,α22,α32,α42,α52,α62,β12,β22,β32,β42,β52,β62,c22,c32,c42,c52,c62 = cache.tab

  if integrator.u_modified
    cache.step = 1
  end

  if cache.step % 2 == 1
    # u1
    @. tmp = dt*fsalfirst
    @. u   = uprev + β11*tmp
    # u2
    f( k,  u, p, t+c21*dt)
    @. tmp = α21*tmp + dt*k
    @. u   = u + β21*tmp
    # u3
    f( k,  u, p, t+c31*dt)
    @. tmp = α31*tmp + dt*k
    @. u   = u + β31*tmp
    # u4
    f( k,  u, p, t+c41*dt)
    @. tmp = α41*tmp + dt*k
    @. u   = u + β41*tmp
    # u5 = u
    f( k,  u, p, t+c51*dt)
    @. tmp = α51*tmp + dt*k
    @. u   = u + β51*tmp

    f( k,  u, p, t+dt)
  else
    # u1
    @. tmp = dt*fsalfirst
    @. u   = uprev + β12*tmp
    # u2
    f( k,  u, p, t+c22*dt)
    @. tmp = α22*tmp + dt*k
    @. u   = u + β22*tmp
    # u3
    f( k,  u, p, t+c32*dt)
    @. tmp = α32*tmp + dt*k
    @. u   = u + β32*tmp
    # u4
    f( k,  u, p, t+c42*dt)
    @. tmp = α42*tmp + dt*k
    @. u   = u + β42*tmp
    # u5 = u
    f( k,  u, p, t+c52*dt)
    @. tmp = α52*tmp + dt*k
    @. u   = u + β52*tmp
    # u6 = u
    f( k,  u, p, t+c62*dt)
    @. tmp = α62*tmp + dt*k
    @. u   = u + β62*tmp

    f( k,  u, p, t+dt)
  end
end

function initialize!(integrator,cache::SSPRK33ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
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
  integrator.destats.nf += 3
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK33Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,stage_limiter!,step_limiter! = cache

  # u1
  @.. u = uprev + dt*fsalfirst
  stage_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
  # u2
  @.. u = (3*uprev + u + dt*k) / 4
  stage_limiter!(u, f, t+dt/2)
  f(k,u,p,t+dt/2)
  # u
  @.. u = (uprev + 2*u + 2*dt*k) / 3
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 3
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK53ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
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
  integrator.destats.nf += 5
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK53Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack α30,α32,α40,α43,α52,α54,β10,β21,β32,β43,β54,c1,c2,c3,c4 = cache.tab

  # u1
  @.. tmp = uprev + β10 * dt * fsalfirst
  stage_limiter!(tmp, f, t+c1*dt)
  f( k,  tmp, p, t+c1*dt)
  # u2 -> stored as u
  @.. u = tmp + β21 * dt * k
  stage_limiter!(u, f, t+c2*dt)
  f( k,  u, p, t+c2*dt)
  # u3
  @.. tmp = α30 * uprev + α32 * u + β32 * dt * k
  stage_limiter!(tmp, f, t+c3*dt)
  f( k,  tmp, p, t+c3*dt)
  # u4
  @.. tmp = α40 * uprev + α43 * tmp + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f( k,  tmp, p, t+c4*dt)
  # u
  @.. u = α52 * u + α54 * tmp + β54 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 5
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK53_2N1ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK53_2N1ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α40,α43,β10,β21,β32,β43,β54,c1,c2,c3,c4 = cache
  #stores in u for all intermediate stages
  # u1
  u = uprev + β10 * dt * integrator.fsalfirst
  k = f(u, p, t+c1*dt)
  # u2
  u = u + β21 * dt * k
  k = f(u, p, t+c2*dt)
  # u3
  u = u + β32 * dt * k
  k = f(u, p, t+c3*dt)
  # u4
  u= α40 * uprev + α43 * u + β43 * dt * k
  k = f(u, p, t+c4*dt)
  # u
  u =  u + β54 * dt * k

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK53_2N1Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK53_2N1Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack α40,α43,β10,β21,β32,β43,β54,c1,c2,c3,c4 = cache.tab
  #stores in u for all intermediate stages
  # u1
  @.. u = uprev + β10 * dt * fsalfirst
  stage_limiter!(u, f, t+c1*dt)
  f( k,  u, p, t+c1*dt)
  # u2
  @.. u = u + β21 * dt * k
  stage_limiter!(u, f, t+c2*dt)
  f( k,  u, p, t+c2*dt)
  # u3
  @.. u = u + β32 * dt * k
  stage_limiter!(u, f, t+c3*dt)
  f( k,  u, p, t+c3*dt)
  # u4
  @.. u = α40 * uprev + α43 * u + β43 * dt * k
  stage_limiter!(u, f, t+c4*dt)
  f( k,  u, p, t+c4*dt)
  # u
  @.. u = u + β54 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 5
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK53_2N2ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK53_2N2ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α30,α32,α50,α54,β10,β21,β32,β43,β54,c1,c2,c3,c4 = cache
  #stores in u for all intermediate stages
  # u1
  u = uprev + β10 * dt * integrator.fsalfirst
  k = f(u, p, t+c1*dt)
  # u2 -> stored as u
  u = u + β21 * dt * k
  k = f(u, p, t+c2*dt)
  # u3
  u = α30 * uprev + α32 * u + β32 * dt * k
  k = f(u, p, t+c3*dt)
  # u4
  u = u + β43 * dt * k
  k = f(u, p, t+c4*dt)
  # u
  u = α50 * uprev + α54 * u + β54 * dt * k

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.destats.nf += 5
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK53_2N2Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK53_2N2Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack α30,α32,α50,α54,β10,β21,β32,β43,β54,c1,c2,c3,c4 = cache.tab

  # u1
  @.. u = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u, f, t+c1*dt)
  f( k,  u, p, t+c1*dt)
  # u2 -> stored as u
  @.. u = u + β21 * dt * k
  stage_limiter!(u, f, t+c2*dt)
  f( k,  u, p, t+c2*dt)
  # u3
  @.. u = α30 * uprev + α32 * u + β32 * dt * k
  stage_limiter!(u, f, t+c3*dt)
  f( k,  u, p, t+c3*dt)
  # u4
  @.. u = u + β43 * dt * k
  stage_limiter!(u, f, t+c4*dt)
  f( k,  u, p, t+c4*dt)
  # u
  @.. u = α50* uprev + α54 * u+ β54 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 5
  f( k,  u, p, t+dt)
end
function initialize!(integrator,cache::SSPRK53_HConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRK53_HConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α30,α32,α40,α41,α43,β10,β21,β32,β43,β54,c1,c2,c3,c4 = cache
  #stores in u for all intermediate stages
  # u1
  tmp = uprev + β10 * dt * integrator.fsalfirst
  k = f(tmp, p, t+c1*dt)
  # u2
  u = tmp + β21 * dt * k
  k = f(u, p, t+c2*dt)
  # u3
  u = α30 * uprev + α32 * u + β32 * dt * k
  k = f(u, p, t+c3*dt)
  # u4
  u = α40 * uprev + α41 * tmp + α43 * u + β43 * dt * k
  k = f(u, p, t+c4*dt)
  # u
  u =  u + β54 * dt * k

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.destats.nf += 5
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::SSPRK53_HCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::SSPRK53_HCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack α30,α32,α40,α41,α43,β10,β21,β32,β43,β54,c1,c2,c3,c4 = cache.tab
  #stores in u for all intermediate stages
  # u1
  @.. tmp = uprev + β10 * dt * fsalfirst
  stage_limiter!(tmp, f, t+c1*dt)
  f( k,  tmp, p, t+c1*dt)
  # u2
  @.. u = tmp + β21 * dt * k
  stage_limiter!(u, f, t+c2*dt)
  f( k,  u, p, t+c2*dt)
  # u3
  @.. u = α30 * uprev + α32 * u + β32 * dt * k
  stage_limiter!(u, f, t+c3*dt)
  f( k,  u, p, t+c3*dt)
  # u4
  @.. u = α40 * uprev + α41 * tmp + α43 * u + β43 * dt * k
  stage_limiter!(u, f, t+c4*dt)
  f( k,  u, p, t+c4*dt)
  # u
  @.. u = u + β54 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 5
  f( k,  u, p, t+dt)
end



function initialize!(integrator,cache::SSPRK63ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
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
  integrator.destats.nf += 6
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK63Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,u₂,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack α40,α41,α43,α62,α65,β10,β21,β32,β43,β54,β65,c1,c2,c3,c4,c5 = cache.tab

  # u1 -> stored as u
  @.. u = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u, f, t+c1*dt)
  f( k,  u, p, t+c1*dt)
  # u2
  @.. u₂ = u + β21 * dt * k
  stage_limiter!(u₂, f, t+c2*dt)
  f(k,u₂,p,t+c2*dt)
  # u3
  @.. tmp = u₂ + β32 * dt * k
  stage_limiter!(tmp, f, t+c3*dt)
  f( k,  tmp, p, t+c3*dt)
  # u4
  @.. tmp = α40 * uprev + α41 * u + α43 * tmp + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f( k,  tmp, p, t+c4*dt)
  # u5
  @.. tmp = tmp + β54 * dt * k
  stage_limiter!(tmp, f, t+c5*dt)
  f( k,  tmp, p, t+c5*dt)
  # u
  @.. u = α62 * u₂ + α65 * tmp + β65 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 6
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK73ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
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
  integrator.destats.nf += 7
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK73Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,u₁,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack α40,α43,α50,α51,α54,α73,α76,β10,β21,β32,β43,β54,β65,β76,c1,c2,c3,c4,c5,c6 = cache.tab

  # u1
  @.. u₁ = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u₁, f, t+c1*dt)
  f(k,u₁,p,t+c1*dt)
  # u2
  @.. tmp = u₁ + β21 * dt * k
  stage_limiter!(tmp, f, t+c2*dt)
  f( k,  tmp, p, t+c2*dt)
  # u3 -> stored as u
  @.. u = tmp + β32 * dt * k
  stage_limiter!(u, f, t+c3*dt)
  f( k,  u, p, t+c3*dt)
  # u4
  @.. tmp = α40 * uprev + α43 * u + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f( k,  tmp, p, t+c4*dt)
  # u5
  @.. tmp = α50 * uprev + α51 * u₁ + α54 * tmp + β54 * dt * k
  stage_limiter!(tmp, f, t+c5*dt)
  f( k,  tmp, p, t+c5*dt)
  # u6
  @.. tmp = tmp + β65 * dt * k
  stage_limiter!(tmp, f, t+c6*dt)
  f( k,  tmp, p, t+c6*dt)
  # u
  @.. u = α73 * u + α76 * tmp + β76 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 7
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK83ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
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
  integrator.destats.nf += 8
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK83Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,u₂,u₃,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack α50,α51,α54,α61,α65,α72,α73,α76,β10,β21,β32,β43,β54,β65,β76,β87,c1,c2,c3,c4,c5,c6,c7 = cache.tab

  # u1 -> save as u
  @.. u = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u, f, t+c1*dt)
  f( k,  u, p, t+c1*dt)
  # u2
  @.. u₂ = u + β21 * dt * k
  stage_limiter!(u₂, f, t+c2*dt)
  f(k,u₂,p,t+c2*dt)
  # u3
  @.. u₃ = u₂ + β32 * dt * k
  stage_limiter!(u₃, f, t+c3*dt)
  f(k,u₃,p,t+c3*dt)
  # u4
  @.. tmp = u₃ + β43 * dt * k
  stage_limiter!(tmp, f, t+c4*dt)
  f( k,  tmp, p, t+c4*dt)
  # u5
  @.. tmp = α50 * uprev + α51 * u + α54 * tmp + β54 * dt * k
  stage_limiter!(tmp, f, t+c5*dt)
  f( k,  tmp, p, t+c5*dt)
  # u6
  @.. tmp = α61 * u + α65 * tmp + β65 * dt * k
  stage_limiter!(tmp, f, t+c6*dt)
  f( k,  tmp, p, t+c6*dt)
  # u7
  @.. tmp = α72 * u₂ + α73 * u₃ + α76 * tmp + β76 * dt * k
  stage_limiter!(tmp, f, t+c7*dt)
  f( k,  tmp, p, t+c7*dt)
  # u
  @.. u = tmp + β87 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 8
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK432ConstantCache)
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

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
  integrator.destats.nf += 4
  if integrator.opts.adaptive
    utilde = utilde - u
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK432Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,utilde,atmp,stage_limiter!,step_limiter! = cache
  dt_2 = dt / 2

  # u1
  @.. u = uprev + dt_2*fsalfirst
  stage_limiter!(u, f, t+dt_2)
  f( k,  u, p, t+dt_2)
  # u2
  @.. u = u + dt_2*k
  stage_limiter!(u, f, t+dt)
  f( k,  u, p, t+dt)
  #
  @.. u = u + dt_2*k
  stage_limiter!(u, f, t+dt+dt_2)
  if integrator.opts.adaptive
    @.. utilde = (uprev + 2*u) / 3
  end
  # u3
  @.. u = (2*uprev + u) / 3
  f( k,  u, p, t+dt_2)
  #
  @.. u = u + dt_2*k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)

  if integrator.opts.adaptive
    @.. utilde = utilde - u
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.destats.nf += 4
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRKMSVS32ConstantCache)
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRKMSVS32ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack u_1,u_2,dts,dtf,μ,v_n = cache

  if integrator.iter == 1
    cache.dts[1] = dt
    cache.dts[2] = dt
    cache.dtf[1] = dt
  end
  accpt = true
  dt = dts[1]


  if cache.step < 3 #starting Procedure
    k = f(u,p,t+dt)
    u = uprev + dt*k
    k = f(u, p, t+dt)
    integrator.destats.nf += 2
    u = (uprev + u + dt*k) / 2
      if cache.step == 1
        u_2 = uprev
      else
        u_1 = uprev
      end
    if integrator.opts.adaptive
      v_n = dt/dts[2]*0.5
      cache.dtf[2] = dtf[1]
      cache.dtf[1] = dt/v_n*0.5
      if v_n > 0.5
        cache.step -= 1
        accpt = false
      end
      cache.dts[3] = dts[2]
      cache.dts[2] = dt
      dt = 0.9*dtf[1]
      μ = min(dtf[1],dtf[2])

    end
  else
    if integrator.opts.adaptive
      Ω = (dts[2] + dts[3])/dt
    else
      Ω = 2
    end
    u = (Ω*Ω - 1)/(Ω*Ω)*(uprev + Ω/(Ω-1)*dt*integrator.fsalfirst) + 1/(Ω*Ω)*u_2
    u_2 = u_1
    u_1 = uprev
    if integrator.opts.adaptive
      v_n = (dts[2]+dts[3]-dt)/(dts[2]+dts[3])*0.5
      dt = (dts[2] + dts[3])/(dts[2]+dts[3]+μ)*μ
      cache.dtf[2] = dtf[1]
      dtf[1]  = dt/v_n*0.5
      cache.dts[3] = dts[2]
      cache.dts[2] = dt
      μ = min(dtf[1],dtf[2])
    end
  end
  if accpt == true
    integrator.fsallast = f(u, p, t+dt)
    integrator.destats.nf += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.u = u
  else
    integrator.fsallast = f(uprev, p, t+dt)
    integrator.destats.nf += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.u = uprev
  end
  cache.dts[1] = dt
  cache.step += 1
  cache.u_1 = u_1
  cache.u_2 = u_2
  cache.μ = μ
end

function initialize!(integrator,cache::SSPRKMSVS32Cache)
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRKMSVS32Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,u_1,u_2,stage_limiter!,step_limiter! = cache

  if cache.step < 3
    @.. u = uprev + dt*fsalfirst
    stage_limiter!(u, f, t+dt)
    f(k,u, p, t+dt)
    integrator.destats.nf += 1
    @.. u = (uprev + u + dt*k) / 2
    stage_limiter!(u, f, t+dt)
    step_limiter!(u, f, t+dt)

    if cache.step == 1
      cache.u_2 .= uprev
    else
      cache.u_1 .= uprev
    end
  else
    Ω = 2
    @.. u = ((Ω*Ω - 1)/(Ω*Ω))*(uprev + (Ω/(Ω-1))*dt*fsalfirst) + (1/(Ω*Ω))*cache.u_2
    cache.u_2 .= u_1
    cache.u_1 .= uprev
    stage_limiter!(u, f, t+dt)
    step_limiter!(u, f, t+dt)
  end
  cache.step += 1
  integrator.destats.nf += 1
  f(k,u, p, t+dt)
end


function initialize!(integrator,cache::SSPRKMSVS43ConstantCache)
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::SSPRKMSVS43ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack u_1,u_2,u_3,k1,k2,k3 = cache

  if cache.step < 4
    u = uprev + dt*integrator.fsalfirst
    k = f(u, p, t+dt)
    u = (uprev + u + dt*k) / 2
      if cache.step == 1
        u_3 = uprev
        cache.k3 = f(u_3,p,t+dt)
        integrator.destats.nf += 1
      end
      if cache.step == 2
        u_2 = uprev
        cache.k2 = f(u_2,p,t+dt)
        integrator.destats.nf += 1
      end
      if cache.step == 3
        u_1 = uprev
        cache.k1 = f(u_1,p,t+dt)
        integrator.destats.nf += 1
      end
  # u
  else
    u = (16/27)*(uprev + 3*dt*integrator.fsalfirst) + (11/27)*(u_3 + (12/11)*dt*k3)
    cache.k3 = k2
    cache.k2 = k1
    cache.k1 = integrator.fsalfirst
    u_3 = u_2
    u_2 = u_1
    u_1 = uprev
  end
  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
  cache.step += 1
  cache.u_1 = u_1
  cache.u_2 = u_2
  cache.u_3 = u_3
end

function initialize!(integrator,cache::SSPRKMSVS43Cache)
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end


@muladd function perform_step!(integrator,cache::SSPRKMSVS43Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,u_1,u_2,u_3,stage_limiter!,step_limiter!,k1,k2,k3 = cache

  if cache.step < 4
    @.. u = uprev + dt*fsalfirst
    stage_limiter!(u, f, t+dt)
    f(k,u, p, t+dt)
    integrator.destats.nf += 1
    @.. u = (uprev + u + dt*k) / 2
    stage_limiter!(u, f, t+dt)
    step_limiter!(u, f, t+dt)
    if cache.step == 1
      cache.u_3 .= uprev
      f(k3,u_3,p,t+dt)
      integrator.destats.nf += 1
    end
    if cache.step == 2
      cache.u_2 .= uprev
      f(k2,u_2,p,t+dt)
      integrator.destats.nf += 1
    end
    if cache.step == 3
      cache.u_1 .= uprev
      f(k1,u_1,p,t+dt)
      integrator.destats.nf += 1
    end
  # u
  else
    @.. u = (16/27)*(uprev + 3*dt*fsalfirst) + (11/27)*(u_3 + (12/11)*dt*k3)
    stage_limiter!(u, f, t+dt)
    step_limiter!(u, f, t+dt)
    cache.k3 .= k2
    cache.k2 .= k1
    cache.k1 .= fsalfirst
    cache.u_3 .= u_2
    cache.u_2 .= u_1
    cache.u_1 .= uprev
  end
  cache.step += 1
  integrator.destats.nf += 1
  f( k, u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK932ConstantCache)
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

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
  integrator.destats.nf += 5
  # u6
  u = u + dt_6*k
  if integrator.opts.adaptive
    k = f(u, p, t+dt)
    integrator.destats.nf += 1
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
  integrator.destats.nf += 4
  if integrator.opts.adaptive
    utilde = utilde - u
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK932Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,utilde,atmp,stage_limiter!,step_limiter! = cache
  dt_6 = dt / 6
  dt_3 = dt / 3
  dt_2 = dt / 2

  # u1
  @.. u = uprev + dt_6*fsalfirst
  stage_limiter!(u, f, t+dt_6)
  f( k,  u, p, t+dt_6)
  # u2
  @.. u = u + dt_6*k
  stage_limiter!(u, f, t+dt_3)
  f( k,  u, p, t+dt_3)
  # u3
  @.. u = u + dt_6*k
  stage_limiter!(u, f, t+dt_2)
  f( k,  u, p, t+dt_2)
  # u4
  @.. u = u + dt_6*k
  stage_limiter!(u, f, t+2*dt_3)
  f( k,  u, p, t+2*dt_3)
  # u5
  @.. u = u + dt_6*k
  stage_limiter!(u, f, t+5*dt_6)
  f( k,  u, p, t+5*dt_6)
  integrator.destats.nf += 5
  # u6
  @.. u = u + dt_6*k
  if integrator.opts.adaptive
    stage_limiter!(u, f, t+dt)
    f( k,  u, p, t+dt)
    integrator.destats.nf += 1
    @.. utilde = (uprev + 6*u + 6*dt*k) / 7
    stage_limiter!(utilde, f, t+dt)
  end
  # u6*
  @.. u = (3*uprev + dt_2*fsalfirst + 2*u) / 5
  stage_limiter!(u, f, t+dt_6)
  f( k,  u, p, t+dt_2)
  # u7*
  @.. u = u + dt_6*k
  stage_limiter!(u, f, t+2*dt_3)
  f( k,  u, p, t+2*dt_3)
  # u8*
  @.. u = u + dt_6*k
  stage_limiter!(u, f, t+5*dt_6)
  f( k,  u, p, t+5*dt_6)
  # u9*
  @.. u = u + dt_6*k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)

  if integrator.opts.adaptive
    @.. utilde = utilde - u
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.destats.nf += 4
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK54ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
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
  integrator.destats.nf += 5
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK54Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,k₃,u₂,u₃,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack β10,α20,α21,β21,α30,α32,β32,α40,α43,β43,α52,α53,β53,α54,β54,c1,c2,c3,c4 = cache.tab

  # u₁
  @.. u₂ = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u₂, f, t+c1*dt)
  f(k,u₂,p,t+c1*dt)
  # u₂
  @.. u₂ = α20 * uprev + α21 * u₂ + β21 * dt * k
  stage_limiter!(u₂, f, t+c2*dt)
  f(k,u₂,p,t+c2*dt)
  # u₃
  @.. u₃ = α30 * uprev + α32 * u₂ + β32 * dt * k
  stage_limiter!(u₃, f, t+c3*dt)
  f(k₃,u₃,p,t+c3*dt)
  # u₄ -> stored as tmp
  @.. tmp = α40 * uprev + α43 * u₃ + β43 * dt * k₃
  stage_limiter!(tmp, f, t+c4*dt)
  f( k,  tmp, p, t+c4*dt)
  # u
  @.. u = α52 * u₂ + α53 * u₃ + β53 * dt * k₃ + α54 * tmp + β54 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 5
  f( k,  u, p, t+dt)
end

function initialize!(integrator,cache::VTSRKConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::VTSRKConstantCache,repeat_step=false)
  @unpack t,dt,uprev2,uprev,u,f,p = integrator
  @unpack β10,α20,α21,β21,α30,α32,β32,α40,α43,β43,α52,α53,β53,α54,β54,d7,n2,n3,n6,n8,q20,q30,q70,q80,q21,q41,q51,q61,q71,q81,q32,q82,q43,q54,q65,q76,q87,r,c2,c3,c4,c5,c6,c7,c8,c1hat,c2hat,c3hat,c4hat = cache

  if integrator.u_modified
      cache.step = 1
    end
  # u₁
  if cache.step == 1
    u₂ = uprev + β10 * dt * integrator.fsalfirst*0.25
    k = f(u₂,p,t+c1hat*dt*0.25)
    # u₂
    u₂ = α20 * uprev + α21 * u₂ + β21 * dt * k * 0.25
    k = f(u₂,p,t+c2hat*dt*0.25)
    # u₃
    u₃ = α30 * uprev + α32 * u₂ + β32 * dt * k * 0.25
    k₃ = f(u₃,p,t+c3hat*dt*0.25)
    # u₄ -> stored as tmp
    tmp = α40 * uprev + α43 * u₃ + β43 * dt * k₃*0.25
    k = f(tmp, p, t+c4hat*dt*0.25)
    # u
    u = α52 * u₂ + α53 * u₃ + β53 * dt * k₃ + α54 * tmp + β54 * dt * k*0.25

    u₃ = u

    y2 = (1 - q20 - q21)*u + q20*(uprev + dt*0.25*f(uprev,p,t)/r) + q21*(u + dt*0.25*f(u,p,t + 0.25*dt)/r)
    y3 = (1 - q30 - q32)*u + q30*(uprev + dt*0.25*f(uprev,p,t)/r) + q32*(y2 + dt*0.25*f(y2,p,t + 0.25*dt + c2*0.25*dt)/r)
    y4 = (1 - q41 - q43)*u + q41*(u + dt*0.25*f(u,p,t + 0.25*dt)/r) + q43*(y3 + dt*0.25*f(y3,p,t + 0.25*dt + c3*0.25*dt)/r)
    y5 = (1 - q51 - q54)*u + q51*(u + dt*0.25*f(u,p,t + 0.25*dt)/r) + q54*(y4 + dt*0.25*f(y4,p,t + 0.25*dt + c4*0.25*dt)/r)
    y6 = (1 - q61 - q65)*u + q61*(u + dt*0.25*f(u,p,t + 0.25*dt)/r) + q65*(y5 + dt*0.25*f(y5,p,t + 0.25*dt + c5*0.25*dt)/r)
    y7 = d7*uprev + (1 - d7 - q70 - q71 - q76)*u + q70*(uprev + dt*0.25*f(uprev,p,t)/r) + q71*(u + dt*0.25*f(u,p,t + 0.25*dt)/r) + q76*(y6 + dt*0.25*f(y6,p,t + 0.25*dt + c6*0.25*dt)/r)
    y8 = (1 - q80 - q81 - q82 - q87)*u + q80*(uprev + dt*0.25*f(uprev,p,t)/r) +  q81*(u + dt*0.25*f(u,p,t + 0.25*dt)/r) + q82*(y2 + dt*0.25*f(y2,p,t + 0.25*dt + c2*0.25*dt)/r) + q87*(y7 + dt*0.25*f(y7,p,t + 0.25*dt + c7*0.25*dt)/r)

    u = (1 - n2 - n3 - n6 - n8)*u + n2*(y2 + dt*0.25*f(y2,p,t + 0.25*dt + c2*0.25*dt)/r) + n3*(y3 + dt*0.25*f(y3,p,t + 0.25*dt + c3*0.25*dt)/r) + n6*(y6 + dt*0.25*f(y6,p,t + 0.25*dt + c6*0.25*dt)/r) + n8*(y8 + dt*0.25*f(y8,p,t + 0.25*dt + c8*0.25*dt)/r)

    y2 = (1 - q20 - q21)*u + q20*(u₃ + dt*0.5*f(u₃,p,t)/r) + q21*(u + dt*0.5*f(u,p,t + 0.5*dt)/r)
    y3 = (1 - q30 - q32)*u + q30*(u₃ + dt*0.5*f(u₃,p,t)/r) + q32*(y2 + dt*0.5*f(y2,p,t + 0.5*dt + c2*0.5*dt)/r)
    y4 = (1 - q41 - q43)*u + q41*(u + dt*0.5*f(u,p,t + 0.5*dt)/r) + q43*(y3 + dt*0.5*f(y3,p,t + 0.5*dt + c3*0.5*dt)/r)
    y5 = (1 - q51 - q54)*u + q51*(u + dt*0.5*f(u,p,t + 0.5*dt)/r) + q54*(y4 + dt*0.5*f(y4,p,t + 0.5*dt + c4*0.5*dt)/r)
    y6 = (1 - q61 - q65)*u + q61*(u + dt*0.5*f(u,p,t + 0.5*dt)/r) + q65*(y5 + dt*0.5*f(y5,p,t + 0.5*dt + c5*0.5*dt)/r)
    y7 = d7*u₃ + (1 - d7 - q70 - q71 - q76)*u + q70*(u₃ + dt*0.5*f(u₃,p,t)/r) + q71*(u + dt*0.5*f(u,p,t + 0.5*dt)/r) + q76*(y6 + dt*0.5*f(y6,p,t + 0.5*dt + c6*0.5*dt)/r)
    y8 = (1 - q80 - q81 - q82 - q87)*u + q80*(u₃ + dt*0.5*f(u₃,p,t)/r) +  q81*(u + dt*0.5*f(u,p,t + 0.5*dt)/r) + q82*(y2 + dt*0.5*f(y2,p,t + 0.5*dt + c2*0.5*dt)/r) + q87*(y7 + dt*0.5*f(y7,p,t + 0.5*dt + c7*0.5*dt)/r)

    u = (1 - n2 - n3 - n6 - n8)*u + n2*(y2 + dt*0.5*f(y2,p,t + 0.5*dt + c2*0.5*dt)/r) + n3*(y3 + dt*0.5*f(y3,p,t + 0.5*dt + c3*0.5*dt)/r) + n6*(y6 + dt*0.5*f(y6,p,t + 0.5*dt + c6*0.5*dt)/r) + n8*(y8 + dt*0.5*f(y8,p,t + 0.5*dt + c8*0.5*dt)/r)

  else
    y2 = (1 - q20 - q21)*uprev + q20*(uprev2 + dt*f(uprev2,p,t-dt)/r) + q21*(uprev + dt*f(uprev,p,t)/r)
    y3 = (1 - q30 - q32)*uprev + q30*(uprev2 + dt*f(uprev2,p,t-dt)/r) + q32*(y2 + dt*f(y2,p,t + c2*dt)/r)
    y4 = (1 - q41 - q43)*uprev + q41*(uprev + dt*f(uprev,p,t)/r) + q43*(y3 + dt*f(y3,p,t + c3*dt)/r)
    y5 = (1 - q51 - q54)*uprev + q51*(uprev + dt*f(uprev,p,t)/r) + q54*(y4 + dt*f(y4,p,t + c4*dt)/r)
    y6 = (1 - q61 - q65)*uprev + q61*(uprev + dt*f(uprev,p,t)/r) + q65*(y5 + dt*f(y5,p,t + c5*dt)/r)
    y7 = d7*uprev2 + (1 - d7 - q70 - q71 - q76)*uprev + q70*(uprev2 + dt*f(uprev2,p,t-dt)/r) + q71*(uprev + dt*f(uprev,p,t)/r) + q76*(y6 + dt*f(y6,p,t + c6*dt)/r)
    y8 = (1 - q80 - q81 - q82 - q87)*uprev + q80*(uprev2 + dt*f(uprev2,p,t-dt)/r) +  q81*(uprev + dt*f(uprev,p,t)/r) + q82*(y2 + dt*f(y2,p,t + c2*dt)/r) + q87*(y7 + dt*f(y7,p,t + c7*dt)/r)

    u = (1 - n2 - n3 - n6 - n8)*uprev + n2*((y2 + dt*f(y2,p,t + c2*dt)/r)) + n3*((y3 + dt*f(y3,p,t + c3*dt)/r)) + n6*(y6 + dt*f(y6,p,t + c6*dt)/r) + n8*(y8 + dt*f(y8,p,t + c8*dt)/r)
  end

  cache.step += 1
  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::VTSRKCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::VTSRKCache,repeat_step=false)
  @unpack t,dt,uprev2,uprev,u,f,p = integrator
  @unpack k,k₃,u₂,u₃,y1,y2,y3,y4,y5,y6,y7,y8,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  @unpack β10,α20,α21,β21,α30,α32,β32,α40,α43,β43,α52,α53,β53,α54,β54,d7,n2,n3,n6,n8,q20,q30,q70,q80,q21,q41,q51,q61,q71,q81,q32,q82,q43,q54,q65,q76,q87,c2,c3,c4,c5,c6,c7,c8,c1hat,c2hat,c3hat,c4hat = cache

  # u₁
  @.. u₂ = uprev + β10 * dt * integrator.fsalfirst
  stage_limiter!(u₂, f, t+c1hat*dt)
  f(k,u₂,p,t+c1*dt)
  # u₂
  @.. u₂ = α20 * uprev + α21 * u₂ + β21 * dt * k
  stage_limiter!(u₂, f, t+c2hat*dt)
  f(k,u₂,p,t+c2*dt)
  # u₃
  @.. u₃ = α30 * uprev + α32 * u₂ + β32 * dt * k
  stage_limiter!(u₃, f, t+c3hat*dt)
  f(k₃,u₃,p,t+c3hat*dt)
  # u₄ -> stored as tmp
  @.. tmp = α40 * uprev + α43 * u₃ + β43 * dt * k₃
  stage_limiter!(tmp, f, t+c4hat*dt)
  f( k,  tmp, p, t+c4hat*dt)
  # u
  @.. u = α52 * u₂ + α53 * u₃ + β53 * dt * k₃ + α54 * tmp + β54 * dt * k
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 5
  f( k,  u, p, t+dt)
end


function initialize!(integrator,cache::SSPRK104ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
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
  integrator.destats.nf += 10
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::SSPRK104Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,k₄,tmp,fsalfirst,stage_limiter!,step_limiter! = cache
  dt_6 = dt/6
  dt_3 = dt/3
  dt_2 = dt/2

  @.. tmp = uprev + dt_6 * integrator.fsalfirst
  stage_limiter!(tmp, f, t+dt_6)
  f( k,  tmp, p, t+dt_6)
  @.. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt_3)
  f( k,  tmp, p, t+dt_3)
  @.. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt_2)
  f( k,  tmp, p, t+dt_2)
  @.. u = tmp + dt_6 * k
  stage_limiter!(u, f, t+2*dt_3)
  f(k₄,u,p,t+2*dt_3)
  @.. tmp = (3*uprev + 2*u + 2*dt_6 * k₄) / 5
  stage_limiter!(tmp, f, t+dt_3)
  f( k,  tmp, p, t+dt_3)
  @.. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt_2)
  f( k,  tmp, p, t+dt_2)
  @.. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+2*dt_3)
  f( k,  tmp, p, t+2*dt_3)
  @.. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+5*dt_6)
  f( k,  tmp, p, t+5*dt_6)
  @.. tmp = tmp + dt_6 * k
  stage_limiter!(tmp, f, t+dt)
  f( k,  tmp, p, t+dt)
  @.. u = (uprev + 9*(u + dt_6*k₄) + 15*(tmp + dt_6*k)) / 25
  stage_limiter!(u, f, t+dt)
  step_limiter!(u, f, t+dt)
  integrator.destats.nf += 10
  f(k, u, p, t+dt)
end
