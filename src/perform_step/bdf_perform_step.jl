function initialize!(integrator, cache::ABDF2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ABDF2ConstantCache, repeat_step=false)
  @unpack t,f,p = integrator
  @unpack uf,κ,tol,dtₙ₋₁ = cache
  alg = unwrap_alg(integrator, true)
  dtₙ, uₙ, uₙ₋₁, uₙ₋₂ = integrator.dt, integrator.u, integrator.uprev, integrator.uprev2

  if integrator.iter == 1 && !integrator.u_modified
    cache.dtₙ₋₁ = dtₙ
    perform_step!(integrator, cache.eulercache, repeat_step)
    cache.fsalfirstprev = integrator.fsalfirst
    return
  end

  # precalculations
  fₙ₋₁ = integrator.fsalfirst
  κtol = κ*tol
  ρ = dtₙ/dtₙ₋₁
  d = 2/3
  ddt = d*dtₙ
  dtmp = ρ^2/3
  d1 = 1+dtmp
  d2 = -dtmp
  d3 = -(ρ-1)/3

  # calculate W
  W = calc_W!(integrator, cache, ddt, repeat_step)

  zₙ₋₁ = dtₙ*fₙ₋₁
  # initial guess
  if alg.extrapolant == :linear
    z = dtₙ*fₙ₋₁
  else # :constant
    z = zero(uₙ)
  end

  tmp = d1*uₙ₋₁ + d2*uₙ₋₂ + d3*zₙ₋₁
  z,η,iter,fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, d, 1, Val{:newton})
  fail_convergence && return

  uₙ = tmp + d*z
  integrator.fsallast = f(uₙ,p,t+dtₙ)

  if integrator.opts.adaptive
    tmp = integrator.fsallast - (1+dtₙ/dtₙ₋₁)*integrator.fsalfirst + (dtₙ/dtₙ₋₁)*cache.fsalfirstprev
    est = (dtₙ₋₁+dtₙ)/6 * tmp
    atmp = calculate_residuals(est, uₙ₋₁, uₙ, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  ################################### Finalize

  cache.dtₙ₋₁ = dtₙ
  cache.ηold = η
  cache.newton_iters = iter
  if integrator.EEst < one(integrator.EEst)
    cache.fsalfirstprev = integrator.fsalfirst
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = uₙ
  return
end

function initialize!(integrator, cache::ABDF2Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::ABDF2Cache, repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack dz,z,k,b,J,W,tmp,atmp,κ,tol,dtₙ₋₁,zₙ₋₁ = cache
  alg = unwrap_alg(integrator, true)
  uₙ,uₙ₋₁,uₙ₋₂,dtₙ = integrator.u,integrator.uprev,integrator.uprev2,integrator.dt

  if integrator.iter == 1 && !integrator.u_modified
    cache.dtₙ₋₁ = dtₙ
    perform_step!(integrator, cache.eulercache, repeat_step)
    cache.fsalfirstprev .= integrator.fsalfirst
    return
  end

  # precalculations
  fₙ₋₁ = integrator.fsalfirst
  κtol = κ*tol
  ρ = dtₙ/dtₙ₋₁
  d = 2/3
  ddt = d*dtₙ
  dtmp = ρ^2/3
  d1 = 1+dtmp
  d2 = -dtmp
  d3 = -(ρ-1)/3

  new_W = calc_W!(integrator, cache, ddt, repeat_step)

  # initial guess
  @. zₙ₋₁ = dtₙ*fₙ₋₁
  if alg.extrapolant == :linear
    @. z = dtₙ*fₙ₋₁
  else # :constant
    fill!(z, zero(eltype(z)))
  end

  @. tmp = d1*uₙ₋₁ + d2*uₙ₋₂ + d3*zₙ₋₁
  z,η,iter,fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, d, 1, Val{:newton}, new_W)
  fail_convergence && return

  @. uₙ = tmp + d*z

  f(integrator.fsallast, uₙ, p, t+dtₙ)
  if integrator.opts.adaptive
    btilde0 = (dtₙ₋₁+dtₙ)/6
    btilde1 = 1+dtₙ/dtₙ₋₁
    btilde2 = dtₙ/dtₙ₋₁
    @. tmp = btilde0*(integrator.fsallast - btilde1*integrator.fsalfirst + btilde2*cache.fsalfirstprev)
    calculate_residuals!(atmp, tmp, uₙ₋₁, uₙ, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter
  cache.dtₙ₋₁ = dtₙ
  if integrator.EEst < one(integrator.EEst)
    @. cache.fsalfirstprev = integrator.fsalfirst
  end
  return
end

# QNDF1

function initialize!(integrator, cache::QNDF1ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::QNDF1ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uprev2,D,temp_D,D2,R,U,dtₙ₋₁,uf,κ,tol = cache
  cnt = integrator.iter
  if cnt == 1
    perform_step!(integrator, cache.eulercache, repeat_step)
    cache.uprev2 = integrator.uprev
    cache.dtₙ₋₁ = dt
    return
  end
  κ = -0.1850
  γₖ = 1//1
  k = 1
  ρ = dt/dtₙ₋₁
  D[1,1] = uprev - uprev2 # backward diff
  temp_D .= D
  if ρ != 1
    U!(k,cache)
    R!(k,ρ,cache)
    D[1,1] = D[1,1] * (R[1,1] * U[1,1])
    uprev2 = prev_u(uprev,k,t,dt,cache)
  end
  c1 = (1-2κ)/(1-κ)
  c2 = (κ)/(1-κ)
  tmp = c1*uprev + c2*uprev2
  # precalculations
  γ = inv(1-κ)
  γdt = γ*dt
  W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  z = dt*integrator.fsalfirst

  z, η, iter, fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, γ, 1, Val{:newton})
  fail_convergence && return
  u = tmp + γ*z

  if integrator.opts.adaptive
    D2!(u,uprev,k,cache)
    utilde = (κ*γₖ + inv(k+1)) * D2[1,1]
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
    if integrator.EEst > one(integrator.EEst)
      cache.D .= temp_D
      return
    end
  end
  cache.dtₙ₋₁ = dt
  cache.uprev2 = uprev
  cache.ηold = η
  cache.newton_iters = iter
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::QNDF1Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

function perform_step!(integrator,cache::QNDF1Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uprev2,D,temp_D,D2,R,U,tmp,z,k,W,dtₙ₋₁,uf,κ,tol,utilde,atmp = cache
  cnt = integrator.iter
  if cnt == 1
    perform_step!(integrator, cache.eulercache, repeat_step)
    cache.uprev2 .= integrator.uprev
    cache.dtₙ₋₁ = dt
    return
  end
  κ = -0.1850
  γₖ = 1//1
  k = 1
  ρ = dt/dtₙ₋₁
  @. D[1,1] = uprev - uprev2 # backward diff
  temp_D .= D
  if ρ != 1
    U!(k,cache)
    R!(k,ρ,cache)
    @. D[1,1] = D[1,1] * (R[1,1] * U[1,1])
    uprev2 = prev_u(uprev,k,t,dt,cache)
  end
  c1 = (1-2κ)/(1-κ)
  c2 = (κ)/(1-κ)
  @. tmp = c1*uprev + c2*uprev2
  # precalculations
  γ = inv(1-κ)
  γdt = γ*dt
  new_W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  @. z = dt*integrator.fsalfirst

  z, η, iter, fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, γ, 1, Val{:newton}, new_W)
  fail_convergence && return
  @. u = tmp + γ*z

  if integrator.opts.adaptive
    D2!(u,uprev,k,cache)
    @. utilde = (κ*γₖ + inv(k+1)) * D2[1,1]
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
    if integrator.EEst > one(integrator.EEst)
      cache.D .= temp_D
      return
    end
  end
  cache.dtₙ₋₁ = dt
  cache.uprev2 .= uprev
  cache.ηold = η
  cache.newton_iters = iter
  f(integrator.fsallast, u, p, t+dt)
end
