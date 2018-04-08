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
  dtₙ, uₙ, uₙ₋₁, uₙ₋₂ = integrator.dt, integrator.u, integrator.uprev, integrator.uprev2

  if integrator.iter == 1 && !integrator.u_modified
    cache.dtₙ₋₁ = dtₙ
    x = perform_step!(integrator, cache.eulercache, repeat_step)
    cache.fsalfirstprev = integrator.fsalfirst
    return x
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
  uf.t = t
  if typeof(uₙ₋₁) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uₙ₋₁)
    W = I - ddt*J
  else
    J = ForwardDiff.derivative(uf,uₙ₋₁)
    W = 1 - ddt*J
  end

  zₙ₋₁ = dtₙ*fₙ₋₁
  # initial guess
  if integrator.alg.extrapolant == :linear
    z = dtₙ*fₙ₋₁
  else # :constant
    z = zero(uₙ)
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dtₙ
  tmp = d1*uₙ₋₁ + d2*uₙ₋₂ + d3*zₙ₋₁
  uₙ = tmp + d*z
  b = dtₙ*f(uₙ, p, tstep) - z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    uₙ = tmp + d*z
    b = dtₙ*f(uₙ, p, tstep) - z
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z = z + dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  uₙ = tmp + d*z
  integrator.fsallast = f(uₙ,p,tstep)

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
  cache.fsalfirstprev = integrator.fsalfirst

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = uₙ
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
  uₙ,uₙ₋₁,uₙ₋₂,dtₙ = integrator.u,integrator.uprev,integrator.uprev2,integrator.dt

  if integrator.iter == 1 && !integrator.u_modified
    cache.dtₙ₋₁ = dtₙ
    x = perform_step!(integrator, cache.eulercache, repeat_step)
    cache.fsalfirstprev = integrator.fsalfirst
    return x
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

  new_W = calc_W!(integrator, cache, ddt, repeat_step, true)

  # initial guess
  @. zₙ₋₁ = dtₙ*fₙ₋₁
  if integrator.alg.extrapolant == :linear
    @. z = dtₙ*fₙ₋₁
  else # :constant
    fill!(z, zero(eltype(z)))
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dtₙ
  @. tmp = d1*uₙ₋₁ + d2*uₙ₋₂ + d3*zₙ₋₁
  @. uₙ = tmp + d*z
  f(k, uₙ, p, tstep)
  @. b = dtₙ*k - z
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  @. z += dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. uₙ = tmp + d*z
    f(k, uₙ, p, tstep)
    @. b = dtₙ*k - z
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. uₙ = tmp + d*z

  f(integrator.fsallast, uₙ, p, t)
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
  @. cache.fsalfirstprev = integrator.fsalfirst
end
