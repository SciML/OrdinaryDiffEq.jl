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
  k = integrator.fsalfirst
  u = @. uprev + dt*k
  k = f(t+dt,u) # For the interpolation, needs k at the updated point
  integrator.fsallast = k
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

@inline function initialize!(integrator,cache::TrapezoidConstantCache,f=integrator.f)
  cache.uhold[1] = integrator.uprev; cache.C[1] = integrator.uprev
  integrator.fsalfirst = f(integrator.t,integrator.uprev)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline function perform_step!(integrator,cache::TrapezoidConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  @unpack uhold,C,rhs,nl_rhs = cache
  C[1] = first(uhold) + (dt/2)*first(integrator.fsalfirst)
  if integrator.iter > 1 && !integrator.u_modified
    uhold[1] = current_extrapolant(t+dt,integrator)
  end # else uhold is previous value.
  rhs.t = t
  rhs.dt = dt
  rhs.a = dt/2
  nlres = integrator.alg.nlsolve(nl_rhs,uhold)
  uhold[1] = nlres[1]
  k = f(t+dt,uhold[1])
  integrator.fsallast = k
  u = uhold[1]
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end
