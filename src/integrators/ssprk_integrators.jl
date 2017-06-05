@inline function initialize!(integrator,cache::SSPRK22ConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

@inline function perform_step!(integrator,cache::SSPRK22ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  k = integrator.fsalfirst
  tmp = @. uprev + dt*k
  k = f(t+dt,tmp)
  u = @. (uprev + tmp + dt*k) / 2
  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::SSPRK22Cache,f=integrator.f)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@inline function perform_step!(integrator,cache::SSPRK22Cache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k,du,tmp,fsalfirst = cache
  @. tmp = @muladd uprev + dt*integrator.fsalfirst
  f(t+dt,tmp,k)
  @. u = (uprev + tmp + dt*k) / 2
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end


@inline function initialize!(integrator,cache::SSPRK33ConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

@inline function perform_step!(integrator,cache::SSPRK33ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  k = integrator.fsalfirst
  tmp = uprev + dt*k
  k = f(t+dt,tmp)
  tmp = (3*uprev + tmp + dt*k) / 4
  k = f(t+dt/2,tmp)
  u = (uprev + 2*tmp + 2*dt*k) / 3
  integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::SSPRK33Cache,f=integrator.f)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

@inline function perform_step!(integrator,cache::SSPRK33Cache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k,du,tmp,fsalfirst = cache
  @. tmp = @muladd uprev + dt*integrator.fsalfirst
  f(t+dt,tmp,k)
  @. tmp = (3*uprev + tmp + dt*k) / 4
  f(t+dt/2,tmp,k)
  @. u = (uprev + 2*tmp + 2*dt*k) / 3
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end


@inline function initialize!(integrator,cache::SSPRK104ConstantCache,f=integrator.f)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
end

@inline function perform_step!(integrator,cache::SSPRK104ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  k = integrator.fsalfirst
  dt_6 = dt/6
  dt_3 = dt/3
  dt_2 = dt/2

  tmp = uprev + dt_6 * k # u₁
  k = f(t+dt_6, tmp)
  tmp = tmp + dt_6 * k # u₂
  k = f(t+dt_3, tmp)
  tmp = tmp + dt_6 * k # u₃
  k = f(t+dt_2, tmp)
  u₄ = tmp + dt_6 * k # u₄
  k₄ = f(t+2*dt_3, u₄)
  tmp = (3*uprev + 2*u₄ + 2*dt_6 * k₄) / 5 # u₅
  k = f(t+dt_3, tmp)
  tmp = tmp + dt_6 * k # u₆
  k = f(t+dt_2, tmp)
  tmp = tmp + dt_6 * k # u₇
  k = f(t+2*dt_3, tmp)
  tmp = tmp + dt_6 * k # u₈
  k = f(t+5*dt_6, tmp)
  tmp = tmp + dt_6 * k # u₉
  k = f(t+dt, tmp)
  u = (uprev + 9*(u₄ + dt_6*k₄) + 15*(tmp + dt_6*k)) / 25

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

@inline function perform_step!(integrator,cache::SSPRK104Cache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack u₄,k,k₄,du,tmp,fsalfirst = cache
  dt_6 = dt/6
  dt_3 = dt/3
  dt_2 = dt/2

  @. tmp = @muladd uprev + dt_6 * integrator.fsalfirst
  f(t+dt_6, tmp, k)
  @. tmp = @muladd tmp + dt_6 * k
  f(t+dt_3, tmp, k)
  @. tmp = @muladd tmp + dt_6 * k
  f(t+dt_2, tmp, k)
  @. u₄ = @muladd tmp + dt_6 * k
  f(t+2*dt_3, u₄, k₄)
  @. tmp = @muladd (3*uprev + 2*u₄ + 2*dt_6 * k₄) / 5
  f(t+dt_3, tmp, k)
  @. tmp = @muladd tmp + dt_6 * k
  f(t+dt_2, tmp, k)
  @. tmp = @muladd tmp + dt_6 * k
  f(t+2*dt_3, tmp, k)
  @. tmp = @muladd tmp + dt_6 * k
  f(t+5*dt_6, tmp, k)
  @. tmp = @muladd tmp + dt_6 * k
  f(t+dt, tmp, k)
  @. u = @muladd (uprev + 9*(u₄ + dt_6*k₄) + 15*(tmp + dt_6*k)) / 25
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end
