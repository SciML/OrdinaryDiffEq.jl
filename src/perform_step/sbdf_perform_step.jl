# SBDF2

function initialize!(integrator, cache::SBDF2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::SBDF2ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p,alg = integrator
  @unpack uprev2,k₁,k₂ = cache
  cnt = integrator.iter
  f1 = integrator.f.f1
  f2 = integrator.f.f2
  du₁ = f1(uprev,p,t)
  du₂ = integrator.fsalfirst - du₁
  if cnt == 1
    tmp = uprev + dt*du₂
  else
    tmp = (4*uprev - uprev2)/3 + (dt/3)*(4*du₂ - 2*k₂)
  end
  # Implicit part
  # precalculations
  γ = 1//1
  if cnt != 1
   γ = 2//3
  end
  γdt = γ*dt
  W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  zprev = dt*du₁
  z = zprev # Constant extrapolation

  nlcache = nlsolve_cache(alg, cache, z, tmp, W, γ, 1, true)
  z,η,iter,fail_convergence = diffeq_nlsolve!(integrator, nlcache, cache, alg.nonlinsolve)
  fail_convergence && return
  u = tmp + γ*z

  cache.uprev2 = uprev
  cache.k₁ = du₁
  cache.k₂ = du₂
  cache.ηold = η
  cache.newton_iters = iter
  integrator.fsallast = f1(u, p, t+dt) + f2(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::SBDF2Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
end

function perform_step!(integrator, cache::SBDF2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,alg = integrator
  @unpack tmp,uprev2,k,k₁,k₂,du₁,z = cache
  cnt = integrator.iter
  f1 = integrator.f.f1
  f2 = integrator.f.f2
  f1(du₁, uprev, p, t)
  # Explicit part
  if cnt == 1
    @. tmp = uprev + dt * (integrator.fsalfirst - du₁)
  else
    @. tmp = (4*uprev - uprev2)/3 + (dt/3)*(4*(integrator.fsalfirst - du₁) - 2*k₂)
  end
  # Implicit part
  # precalculations
  γ = 1//1
  if cnt != 1
   γ = 2//3
  end
  γdt = γ*dt
  new_W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  @. z = dt*du₁

  nlcache = nlsolve_cache(alg, cache, z, tmp, γ, 1, new_W)
  z,η,iter,fail_convergence = diffeq_nlsolve!(integrator, nlcache, cache, alg.nonlinsolve)
  fail_convergence && return
  @. u = tmp + γ*z

  cache.uprev2 .= uprev
  cache.k₁ .= du₁
  @. cache.k₂ = integrator.fsalfirst - du₁
  cache.ηold = η
  cache.newton_iters = iter
  integrator.f(k,u,p,t+dt)
end

# SBDF3

function initialize!(integrator, cache::SBDF3ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::SBDF3ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p,alg = integrator
  @unpack uprev2,uprev3,k₁,k₂ = cache
  cnt = integrator.iter
  f1 = integrator.f.f1
  f2 = integrator.f.f2
  du₁ = f1(uprev,p,t)
  du₂ = integrator.fsalfirst - du₁
  if cnt == 1
    tmp = uprev + dt*du₂
  elseif cnt == 2
    tmp = (4*uprev - uprev2)/3 + (dt/3)*(4*du₂ - 2*k₁)
  else
    tmp = (6//11) * (3*uprev - 3//2*uprev2 + 1//3*uprev3 + dt*(3*(du₂ - k₁) + k₂))
  end
  # Implicit part
  # precalculations
  if cnt == 1
    γ = 1//1
  elseif cnt == 2
    γ = 2//3
  else
    γ = 6//11
  end
  γdt = γ*dt
  W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  zprev = dt*du₁
  z = zprev # Constant extrapolation

  nlcache = nlsolve_cache(alg, cache, z, tmp, W, γ, 1, true)
  z,η,iter,fail_convergence = diffeq_nlsolve!(integrator, nlcache, cache, alg.nonlinsolve)
  fail_convergence && return
  u = tmp + γ*z

  cache.uprev3 = uprev2
  cache.uprev2 = uprev
  cache.k₂ = k₁
  cache.k₁ = du₂
  cache.ηold = η
  cache.newton_iters = iter
  integrator.fsallast = f1(u, p, t+dt) + f2(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::SBDF3Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
end

function perform_step!(integrator, cache::SBDF3Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,alg = integrator
  @unpack tmp,uprev2,uprev3,k,k₁,k₂,du₁,z = cache
  cnt = integrator.iter
  f1 = integrator.f.f1
  f2 = integrator.f.f2
  f1(du₁, uprev, p, t)
  # Explicit part
  if cnt == 1
    @. tmp = uprev + dt*du₂
  elseif cnt == 2
    @. tmp = (4*uprev - uprev2)/3 + (dt/3)*(4*du₂ - 2*k₁)
  else
    @. tmp = (6//11) * (3*uprev - 3//2*uprev2 + 1//3*uprev3 + dt*(3*((integrator.fsalfirst - du₁) - k₁) + k₂))
  end
  # Implicit part
  # precalculations
  if cnt == 1
    γ = 1//1
  elseif cnt == 2
    γ = 2//3
  else
    γ = 6//11
  end
  γdt = γ*dt
  new_W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  @. z = dt*du₁

  nlcache = nlsolve_cache(alg, cache, z, tmp, γ, 1, new_W)
  z,η,iter,fail_convergence = diffeq_nlsolve!(integrator, nlcache, cache, alg.nonlinsolve)
  fail_convergence && return
  @. u = tmp + γ*z

  cache.uprev3 .= uprev2
  cache.uprev2 .= uprev
  cache.k₂ .= k₁
  @. cache.k₁ = integrator.fsalfirst - du₁
  cache.ηold = η
  cache.newton_iters = iter
  integrator.f(k,u,p,t+dt)
end
