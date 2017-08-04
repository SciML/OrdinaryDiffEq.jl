@inline function initialize!(integrator,cache::SSPRK22ConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@inline @muladd function perform_step!(integrator,cache::SSPRK22ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  tmp = @. uprev + dt*integrator.fsalfirst
  k = f(t+dt,tmp)
  u = @. (uprev + tmp + dt*k) / 2
  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::SSPRK22Cache,f=integrator.f)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@inline @muladd function perform_step!(integrator,cache::SSPRK22Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack k,tmp,fsalfirst = cache
  @. tmp = uprev + dt*integrator.fsalfirst
  f(t+dt,tmp,k)
  @. u = (uprev + tmp + dt*k) / 2
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::SSPRK33ConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@inline @muladd function perform_step!(integrator,cache::SSPRK33ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  tmp = @. uprev + dt*integrator.fsalfirst
  k = f(t+dt,tmp)
  tmp = @. (3*uprev + tmp + dt*k) / 4
  k = f(t+dt/2,tmp)
  u = @. (uprev + 2*tmp + 2*dt*k) / 3
  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::SSPRK33Cache,f=integrator.f)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@inline @muladd function perform_step!(integrator,cache::SSPRK33Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack k,tmp,fsalfirst = cache
  @. tmp = uprev + dt*fsalfirst
  f(t+dt,tmp,k)
  @. tmp = (3*uprev + tmp + dt*k) / 4
  f(t+dt/2,tmp,k)
  @. u = (uprev + 2*tmp + 2*dt*k) / 3
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::SSPRK432ConstantCache,f=integrator.f)
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@inline @muladd function perform_step!(integrator,cache::SSPRK432ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
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
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::SSPRK432Cache,f=integrator.f)
  integrator.kshortsize = 1
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  f(integrator.t, integrator.uprev, integrator.fsalfirst) # Pre-start fsal
end

@inline @muladd function perform_step!(integrator,cache::SSPRK432Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack k,tmp,fsalfirst,utilde,atmp = cache
  dt_2 = dt / 2

   # u1
  @. tmp = uprev + dt_2*fsalfirst
  f(t+dt_2, tmp, k)
  # u2
  @. tmp = tmp + dt_2*k
  f(t+dt, tmp, k)
  #
  @. tmp = tmp + dt_2*k
  if integrator.opts.adaptive
    @. utilde = (uprev + 2*tmp) / 3
  end
  # u3
  @. tmp = (2*uprev + tmp) / 3
  f(t+dt_2, tmp, k)
  #
  @. u = tmp + dt_2*k

  if integrator.opts.adaptive
    @. atmp = (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  f(t+dt, u, k)
  @pack integrator = t,dt,u
end


@inline function initialize!(integrator,cache::SSPRK104ConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline @muladd function perform_step!(integrator,cache::SSPRK104ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
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
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::SSPRK104Cache,f=integrator.f)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@inline @muladd function perform_step!(integrator,cache::SSPRK104Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack k,k₄,tmp,fsalfirst = cache
  dt_6 = dt/6
  dt_3 = dt/3
  dt_2 = dt/2

  @. tmp = uprev + dt_6 * integrator.fsalfirst
  f(t+dt_6, tmp, k)
  @. tmp = tmp + dt_6 * k
  f(t+dt_3, tmp, k)
  @. tmp = tmp + dt_6 * k
  f(t+dt_2, tmp, k)
  @. u = tmp + dt_6 * k
  f(t+2*dt_3, u, k₄)
  @. tmp = (3*uprev + 2*u + 2*dt_6 * k₄) / 5
  f(t+dt_3, tmp, k)
  @. tmp = tmp + dt_6 * k
  f(t+dt_2, tmp, k)
  @. tmp = tmp + dt_6 * k
  f(t+2*dt_3, tmp, k)
  @. tmp = tmp + dt_6 * k
  f(t+5*dt_6, tmp, k)
  @. tmp = tmp + dt_6 * k
  f(t+dt, tmp, k)
  @. u = (uprev + 9*(u + dt_6*k₄) + 15*(tmp + dt_6*k)) / 25
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end
