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
  @unpack t,f,p,uprev,uprev2,dt = integrator
  @unpack uf,κ,tol,dtₙ₊₁ = cache
  uₙ₊₁,uₙ,dtₙ₊₂ = uprev,uprev2,dt

  if integrator.iter == 1 && !integrator.u_modified
    cache.dtₙ₊₁ = dtₙ₊₂
    return perform_step!(integrator, cache.eulercache, repeat_step)
  end

  # precalculations
  κtol = κ*tol
  w = dtₙ₊₂/dtₙ₊₁
  dtmp = (1+2w)
  d = (1+w)/dtmp
  d1 = (1+w)^2/dtmp
  d2 = -w^2/dtmp
  ddt = d*dtₙ₊₂

  # calculate W
  uf.t = t
  if typeof(uₙ₊₁) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uₙ₊₁)
    W = I - ddt*J
  else
    J = ForwardDiff.derivative(uf,uₙ₊₁)
    W = 1 - ddt*J
  end

  # initial guess
  if integrator.alg.extrapolant == :linear
    z = dtₙ₊₂*integrator.fsalfirst
  else # :constant
    z = zero(uₙ₊₂)
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dtₙ₊₂
  tmp = d1*uₙ₊₁ + d2*uₙ
  uₙ₊₂ = tmp + d*z
  b = dtₙ₊₂*f(uₙ₊₂, p, tstep) - z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    uₙ₊₂ = tmp + d*z
    b = dtₙ₊₂*f(uₙ₊₂, p, tstep) - z
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

  uₙ₊₂ = tmp + d*z
  integrator.fsallast = z./dtₙ₊₂

  if integrator.opts.adaptive
    tmp = integrator.fsallast - (1+dtₙ₊₂/dtₙ₊₁)*integrator.fsalfirst + (dtₙ₊₂/dtₙ₊₁)*integrator.k[1]
    est = (dtₙ₊₁+dtₙ₊₂)/6 * tmp
    atmp = calculate_residuals(est, uₙ₊₁, uₙ₊₂, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  ################################### Finalize

  cache.dtₙ₊₁ = dtₙ₊₂
  cache.ηold = η
  cache.newton_iters = iter

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = uₙ₊₂
end
