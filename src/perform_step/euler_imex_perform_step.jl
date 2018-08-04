function initialize!(integrator, cache::IMEXEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::IMEXEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,nlsolve = cache
  nlsolve!, nlcache = nlsolve, nlsolve.cache
  f1 = integrator.f.f1
  du₁ = f1(uprev,p,t)
  du₂ = integrator.fsalfirst - du₁
  # Explicit Part
  nlcache.tmp = uprev + dt*du₂
  # Implicit part
  # Precalculations
  typeof(nlsolve!) <: NLNewton && ( nlcache.W = calc_W!(integrator, cache, dt, repeat_step) )

  # initial guess
  alg = unwrap_alg(integrator, true)
  if alg.extrapolant == :linear
    z = dt*du₁
  else # :constant
    z = zero(u)
  end
  nlcache.z = z

  z,η,iter,fail_convergence = nlsolve!(integrator)
  fail_convergence && return
  u = nlcache.tmp + z

  nlcache.ηold = η
  nlcache.nl_iters = iter
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::IMEXEulerCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

function perform_step!(integrator, cache::IMEXEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack du₁,z,W,tmp,nlsolve = cache
  nlsolve!, nlcache = nlsolve, nlsolve.cache
  f1 = integrator.f.f1
  f1(du₁,uprev,p,t)
  # Explicit part
  @. tmp = uprev + dt*(integrator.fsalfirst - du₁)
  # Implicit part
  # Precalculations
  typeof(nlsolve!) <: NLNewton && ( calc_W!(integrator, cache, dt, repeat_step) )

  # initial guess
  alg = unwrap_alg(integrator, true)
  if alg.extrapolant == :linear
    @. z = dt*integrator.fsalfirst
  else # :constant
    z .= zero(eltype(u))
  end

  z,η,iter,fail_convergence = nlsolve!(integrator)
  fail_convergence && return
  @. u = tmp + z

  nlcache.ηold = η
  nlcache.nl_iters = iter
  f(integrator.fsallast,u,p,t+dt)
end
