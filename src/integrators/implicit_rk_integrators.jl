@inline function initialize!(integrator,cache::ImplicitEulerConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline function perform_step!(integrator,cache::ImplicitEulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  @unpack uf = cache
  uf.t = t
  if integrator.iter > 1 && !integrator.u_modified
    u = current_extrapolant(t+dt,integrator)
  else
    u = uprev
  end
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*J
  end
  z = u - uprev
  iter = 0
  do_newton = true
  κ = 1e-2
  tol = 1e-8

  iter += 1
  b = -z + dt*f(t+dt,uprev + z)
  dz = J\b
  ndz = abs(dz)

  while do_newton
    iter += 1
    b = -z + dt*f(t+dt,uprev + z)
    dz = J\b
    ndzprev = ndz
    ndz = abs(dz)
    θ = ndz/ndzprev
    η = θ/(1-θ)
    if η*ndz <= κ*tol
      do_newton=false
    end
    z = z + dz
  end

  u = uprev + z
  integrator.fsallast = f(t+dt,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::ImplicitEulerCache,f=integrator.f)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline function perform_step!(integrator,cache::ImplicitEulerCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack C,dual_cache,k,nl_rhs,rhs,uhold = cache
  copy!(C,uhold)
  if integrator.iter > 1 && !integrator.u_modified
    current_extrapolant!(u,t+dt,integrator)
  end # else uhold is previous value.
  rhs.t = t
  rhs.dt = dt
  rhs.a = dt
  nlres = integrator.alg.nlsolve(nl_rhs,uhold)
  copy!(uhold,nlres)
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::TrapezoidConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline function perform_step!(integrator,cache::TrapezoidConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  @unpack uf = cache
  uf.t = t
  dto2 = dt/2

  if integrator.iter > 1 && !integrator.u_modified
    u = current_extrapolant(t+dt,integrator)
  else
    u = uprev
  end
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dto2*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dto2*J
  end
  z = u - uprev
  iter = 0
  do_newton = true
  κ = 1e-2
  tol = 1e-8

  iter += 1
  b = -z + dto2*f(t+dt,uprev + z)
  dz = J\b
  ndz = abs(dz)

  while do_newton
    iter += 1
    b = -z + dto2*f(t+dto2,uprev + z)
    dz = J\b
    ndzprev = ndz
    ndz = abs(dz)
    θ = ndz/ndzprev
    η = θ/(1-θ)
    if η*ndz <= κ*tol
      do_newton=false
    end
    z = z + dz
  end

  u = uprev + 2z
  integrator.fsallast = f(t+dt,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::TrapezoidCache,f=integrator.f)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = cache.k
  f(integrator.t,integrator.uprev,integrator.fsalfirst)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline function perform_step!(integrator,cache::TrapezoidCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack C,dual_cache,k,rhs,nl_rhs,uhold = cache
  C .= uhold .+ (dt/2).*vec(integrator.fsalfirst)
  if integrator.iter > 1 && !integrator.u_modified
    current_extrapolant!(u,t+dt,integrator)
  end # else uhold is previous value.
  # copy!(rhs.fsalfirst,fsalfirst) Implicitly done by pointers: fsalfirst === fsalfirst == rhs.fsalfirst
  rhs.t = t
  rhs.dt = dt
  rhs.a = dt/2
  nlres = integrator.alg.nlsolve(nl_rhs,uhold)
  copy!(uhold,nlres)
  f(t+dt,u,k)
  @pack integrator = t,dt,u
end
