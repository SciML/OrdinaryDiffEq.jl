function initialize!(integrator, cache::RadauIIA5ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::RadauIIA5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack T11, T12, T13, T21, T22, T23, T31, TI11, TI12, TI13, TI21, TI22, TI23, TI31, TI32, TI33 = cache.tab
  @unpack c1, c2, γ, α, β, e1, e2, e3 = cache.tab
  @unpack tol, κ = cache
  @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
  alg = unwrap_alg(integrator, true)
  @unpack min_iter, max_iter = alg
  mass_matrix = integrator.f.mass_matrix
  is_compos = integrator.alg isa CompositeAlgorithm

  # precalculations
  κtol = κ*tol # used in Newton iteration
  γdt, αdt, βdt = γ/dt, α/dt, β/dt
  J = calc_J(integrator, cache, is_compos)
  if u isa Number
    LU1 = γdt*mass_matrix - J
    LU2 = (αdt + βdt*im)*mass_matrix - J
  else
    LU1 = lu(γdt*mass_matrix - J)
    LU2 = lu((αdt + βdt*im)*mass_matrix - J)
  end

  # TODO better initial guess
  z1 = w1 = map(zero, u)
  z2 = w2 = map(zero, u)
  z3 = w3 = map(zero, u)

  # Newton iteration
  local nw, η
  do_newton = true
  iter = 1
  while do_newton && iter < max_iter
    ff1 = f(uprev+z1, p, t+c1*dt)
    ff2 = f(uprev+z2, p, t+c2*dt)
    ff3 = f(uprev+z3, p, t+   dt) # c3 = 1

    # transform `u'` to `z`
    z1 = @. TI11 * ff1 + TI12 * ff2 + TI13 * ff3
    z2 = @. TI21 * ff1 + TI22 * ff2 + TI23 * ff3
    z3 = @. TI31 * ff1 + TI32 * ff2 + TI33 * ff3

    if mass_matrix isa UniformScaling # `UniformScaling` doesn't play nicely with broadcast
      Mw1 = @. mass_matrix.λ * w1
      Mw2 = @. mass_matrix.λ * w2
      Mw3 = @. mass_matrix.λ * w3
    else
      Mw1 = mass_matrix*w1
      Mw2 = mass_matrix*w2
      Mw3 = mass_matrix*w3
    end

    z1 = LU1 \ (@. z1 - γdt*Mw1)
    z2 = @. z2 - αdt*Mw2 + βdt*Mw3
    z3 = @. z3 - βdt*Mw2 - αdt*Mw3
    z23 = LU2 \ (@. z2 + z3*im)
    z2 = real(z23)
    z3 = imag(z23)

    iter != 1 && (nwprev = nw)
    nw = internalnorm(z1) + internalnorm(z2) + internalnorm(z3)

    # check early stopping criterion
    if iter != 1
      θ = nw/nwprev
      θ ≥ 1 || nw * θ^(max_iter - iter) > κtol * (1 - θ) && break
    end

    w1 = @. w1 + z1
    w2 = @. w2 + z2
    w3 = @. w3 + z3

    # check stopping criterion
    if iter == 1
      η = max(cache.ηold, eps(eltype(reltol)))^(0.8)
      do_newton = iszero(integrator.success_iter) || iter < min_iter || η * nw > κtol
    else
      η = θ / (1 - θ) # calculated for possible early stopping
      do_newton = iter < min_iter || η * nw > κtol
    end

    # transform `w` to `z`
    z1 = @. T11 * w1 + T12 * w2 + T13 * w3
    z2 = @. T21 * w1 + T22 * w2 + T23 * w3
    z3 = @. T31 * w1 +       w2           # T32 = 1, T33 = 0
    iter += 1
  end
  cache.ηold = η
  integrator.force_stepfail = do_newton

  u = @. uprev + z3

  if adaptive
    #tmp = r*integrator.opts.internalnorm.((u - uprev)/dt1 - (uprev - uprev2)/dt2)
    atmp = calculate_residuals(tmp, uprev, u, abstol, reltol, internalnorm)
    integrator.EEst = internalnorm(atmp)
  end

  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end
