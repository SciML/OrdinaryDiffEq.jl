function initialize!(integrator, cache::RadauIIA5ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  nothing
end

function initialize!(integrator, cache::RadauIIA5Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
  if integrator.opts.adaptive
    @unpack abstol, reltol = integrator.opts
    if reltol isa Number
      cache.rtol = reltol^(2/3) / 10
      cache.atol = cache.rtol * (abstol / reltol)
    else
      @. cache.rtol = reltol^(2/3) / 10
      @. cache.atol = cache.rtol * (abstol / reltol)
    end
  end
  nothing
end

@muladd function perform_step!(integrator, cache::RadauIIA5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack T11, T12, T13, T21, T22, T23, T31, TI11, TI12, TI13, TI21, TI22, TI23, TI31, TI32, TI33 = cache.tab
  @unpack c1, c2, γ, α, β, e1, e2, e3 = cache.tab
  @unpack tol, κ, cont1, cont2, cont3 = cache
  @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
  alg = unwrap_alg(integrator, true)
  @unpack min_iter, max_iter = alg
  mass_matrix = integrator.f.mass_matrix
  is_compos = integrator.alg isa CompositeAlgorithm

  # precalculations
  c1m1 = c1-1
  c2m1 = c2-1
  c1mc2= c1-c2
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
  integrator.destats.nw += 1

  # TODO better initial guess
  if integrator.iter == 1 || integrator.u_modified || alg.extrapolant == :constant
    cache.dtprev = one(cache.dtprev)
    z1 = w1 = map(zero, u)
    z2 = w2 = map(zero, u)
    z3 = w3 = map(zero, u)
    cache.cont1 = map(zero, u)
    cache.cont2 = map(zero, u)
    cache.cont3 = map(zero, u)
  else
    c3′ = dt/cache.dtprev
    c1′ = c1*c3′
    c2′ = c2*c3′
    z1  = @. c1′ * (cont1 + (c1′-c2m1) * (cont2 + (c1′-c1m1) * cont3))
    z2  = @. c2′ * (cont1 + (c2′-c2m1) * (cont2 + (c2′-c1m1) * cont3))
    z3  = @. c3′ * (cont1 + (c3′-c2m1) * (cont2 + (c3′-c1m1) * cont3))
    w1  = @. TI11 * z1 + TI12 * z2 + TI13 * z3
    w2  = @. TI21 * z1 + TI22 * z2 + TI23 * z3
    w3  = @. TI31 * z1 + TI32 * z2 + TI33 * z3
  end

  # Newton iteration
  local ndw, η
  do_newton = true
  iter = 0
  while do_newton && iter < max_iter
    iter += 1
    ff1 = f(uprev+z1, p, t+c1*dt)
    ff2 = f(uprev+z2, p, t+c2*dt)
    ff3 = f(uprev+z3, p, t+   dt) # c3 = 1
    integrator.destats.nf += 3

    fw1 = @. TI11 * ff1 + TI12 * ff2 + TI13 * ff3
    fw2 = @. TI21 * ff1 + TI22 * ff2 + TI23 * ff3
    fw3 = @. TI31 * ff1 + TI32 * ff2 + TI33 * ff3

    if mass_matrix isa UniformScaling # `UniformScaling` doesn't play nicely with broadcast
      Mw1 = @. mass_matrix.λ * w1
      Mw2 = @. mass_matrix.λ * w2
      Mw3 = @. mass_matrix.λ * w3
    else
      Mw1 = mass_matrix*w1
      Mw2 = mass_matrix*w2
      Mw3 = mass_matrix*w3
    end

    rhs1 = @. fw1 - γdt*Mw1
    rhs2 = @. fw2 - αdt*Mw2 + βdt*Mw3
    rhs3 = @. fw3 - βdt*Mw2 - αdt*Mw3
    dw1 = LU1 \ rhs1
    dw23 = LU2 \ (@. rhs2 + rhs3*im)
    integrator.destats.nsolve += 2
    dw2 = real(dw23)
    dw3 = imag(dw23)

    iter != 1 && (ndwprev = ndw)
    ndw = internalnorm(dw1,t) + internalnorm(dw2,t) + internalnorm(dw3,t)

    # check early stopping criterion
    if iter != 1
      θ = ndw/ndwprev
      if θ ≥ 1 || ndw * θ^(max_iter - iter) > κtol * (1 - θ)
        break
      end
    end

    w1 = @. w1 + dw1
    w2 = @. w2 + dw2
    w3 = @. w3 + dw3

    # transform `w` to `z`
    z1 = @. T11 * w1 + T12 * w2 + T13 * w3
    z2 = @. T21 * w1 + T22 * w2 + T23 * w3
    z3 = @. T31 * w1 +       w2           # T32 = 1, T33 = 0

    # check stopping criterion
    if iter == 1
      η = max(cache.ηold, eps(eltype(reltol)))^(0.8)
      do_newton = iszero(integrator.success_iter) || iter < min_iter || η * ndw > κtol
    else
      η = θ / (1 - θ) # calculated for possible early stopping
      do_newton = iter < min_iter || η * ndw > κtol
    end
  end
  cache.ηold = η
  integrator.force_stepfail = do_newton
  cache.nl_iters = iter
  do_newton && return

  u = @. uprev + z3

  if adaptive
    e1dt, e2dt, e3dt = e1/dt, e2/dt, e3/dt
    tmp = @. e1dt*z1 + e2dt*z2 + e3dt*z3
    mass_matrix != I && (tmp = mass_matrix*tmp)
    utilde = @. integrator.fsalfirst + tmp
    alg.smooth_est && (utilde = LU1 \ utilde; integrator.destats.nsolve += 1)
    # RadauIIA5 needs a transformed rtol and atol see
    # https://github.com/luchr/ODEInterface.jl/blob/0bd134a5a358c4bc13e0fb6a90e27e4ee79e0115/src/radau5.f#L399-L421
    rtol = @. reltol^(2/3) / 10
    atol = @. rtol * (abstol / reltol)
    atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm,t)
    integrator.EEst = internalnorm(atmp,t)

    if !(integrator.EEst < oneunit(integrator.EEst)) && integrator.iter == 1 || integrator.u_modified
      f0 = f(uprev .+ utilde, p, t)
      integrator.destats.nf += 1
      utilde = @. f0 + tmp
      alg.smooth_est && (utilde = LU1 \ utilde; integrator.destats.nsolve += 1)
      atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm,t)
      integrator.EEst = internalnorm(atmp,t)
    end
  end

  if integrator.EEst <= oneunit(integrator.EEst)
    cache.dtprev = dt
    if alg.extrapolant != :constant
      cache.cont1  = @. (z2 - z3)/c2m1
      tmp          = @. (z1 - z2)/c1mc2
      cache.cont2  = @. (tmp - cache.cont1) / c1m1
      cache.cont3  = @. cache.cont2 - (tmp - z1/c1)/c2
    end
  end

  integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
  return
end

@muladd function perform_step!(integrator, cache::RadauIIA5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,fsallast,fsalfirst = integrator
  @unpack T11, T12, T13, T21, T22, T23, T31, TI11, TI12, TI13, TI21, TI22, TI23, TI31, TI32, TI33 = cache.tab
  @unpack c1, c2, γ, α, β, e1, e2, e3 = cache.tab
  @unpack tol, κ, cont1, cont2, cont3 = cache
  @unpack z1, z2, z3, w1, w2, w3,
          dw1, dw23,
          k, k2, k3, fw1, fw2, fw3,
          J, W1, W2,
          tmp, atmp, jac_config, linsolve1, linsolve2, rtol, atol = cache
  @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
  alg = unwrap_alg(integrator, true)
  @unpack min_iter, max_iter = alg
  mass_matrix = integrator.f.mass_matrix
  is_compos = integrator.alg isa CompositeAlgorithm

  # precalculations
  c1m1 = c1-1
  c2m1 = c2-1
  c1mc2= c1-c2
  κtol = κ*tol # used in Newton iteration
  γdt, αdt, βdt = γ/dt, α/dt, β/dt
  new_W = true
  if repeat_step || (alg_can_repeat_jac(alg) &&
                     (!integrator.last_stepfail && cache.nl_iters == 1 &&
                      cache.ηold < alg.new_jac_conv_bound))
    new_jac = false
  else
    new_jac = true
    calc_J!(integrator, cache, is_compos)
  end
  # skip calculation of W if step is repeated
  if !repeat_step && (!alg_can_repeat_jac(alg) ||
                      (integrator.iter < 1 || new_jac ||
                       abs(dt - (t-integrator.tprev)) > 100eps(typeof(integrator.t))))
    @inbounds for II in CartesianIndices(J)
      W1[II] = γdt * mass_matrix[Tuple(II)...] - J[II]
      W2[II] = (αdt + βdt*im) * mass_matrix[Tuple(II)...] - J[II]
    end
  else
    new_W = false
  end
  new_W && (integrator.destats.nw += 1)

  # TODO better initial guess
  if integrator.iter == 1 || integrator.u_modified || alg.extrapolant == :constant
    cache.dtprev = one(cache.dtprev)
    uzero = zero(eltype(u))
    @. z1 = uzero
    @. z2 = uzero
    @. z3 = uzero
    @. w1 = uzero
    @. w2 = uzero
    @. w3 = uzero
    @. cache.cont1 = uzero
    @. cache.cont2 = uzero
    @. cache.cont3 = uzero
  else
    c3′ = dt/cache.dtprev
    c1′ = c1*c3′
    c2′ = c2*c3′
    @. z1  = c1′ * (cont1 + (c1′-c2m1) * (cont2 + (c1′-c1m1) * cont3))
    @. z2  = c2′ * (cont1 + (c2′-c2m1) * (cont2 + (c2′-c1m1) * cont3))
    @. z3  = c3′ * (cont1 + (c3′-c2m1) * (cont2 + (c3′-c1m1) * cont3))
    @. w1  = TI11 * z1 + TI12 * z2 + TI13 * z3
    @. w2  = TI21 * z1 + TI22 * z2 + TI23 * z3
    @. w3  = TI31 * z1 + TI32 * z2 + TI33 * z3
  end

  # Newton iteration
  local ndw, η
  do_newton = true
  iter = 0
  while do_newton && iter < max_iter
    iter += 1
    @. tmp = uprev + z1
    f(fsallast, tmp, p, t+c1*dt)
    @. tmp = uprev + z2
    f(k2, tmp, p, t+c2*dt)
    @. tmp = uprev + z3
    f(k3, tmp, p, t+   dt) # c3 = 1
    integrator.destats.nf += 3

    @. fw1 = TI11 * fsallast + TI12 * k2 + TI13 * k3
    @. fw2 = TI21 * fsallast + TI22 * k2 + TI23 * k3
    @. fw3 = TI31 * fsallast + TI32 * k2 + TI33 * k3

    if mass_matrix == I
      Mw1 = w1
      Mw2 = w2
      Mw3 = w3
    elseif mass_matrix isa UniformScaling
      mul!(z1, mass_matrix.λ, w1)
      mul!(z2, mass_matrix.λ, w2)
      mul!(z3, mass_matrix.λ, w3)
      Mw1 = z1
      Mw2 = z2
      Mw3 = z3
    else
      mul!(z1, mass_matrix, w1)
      mul!(z2, mass_matrix, w2)
      mul!(z3, mass_matrix, w3)
      Mw1 = z1
      Mw2 = z2
      Mw3 = z3
    end

    @. dw1 = fw1 - γdt*Mw1
    needfactor = iter==1 && new_W
    linsolve1(vec(dw1), W1, vec(dw1), needfactor)
    @. dw23 = complex(fw2 - αdt*Mw2 + βdt*Mw3, fw3 - βdt*Mw2 - αdt*Mw3)
    linsolve2(vec(dw23), W2, vec(dw23), needfactor)
    integrator.destats.nsolve += 2
    dw2 = z2; dw3 = z3
    @. dw2 = real(dw23)
    @. dw3 = imag(dw23)

    iter != 1 && (ndwprev = ndw)
    ndw = internalnorm(dw1,t) + internalnorm(dw2,t) + internalnorm(dw3,t)

    # check early stopping criterion
    if iter != 1
      θ = ndw/ndwprev
      if θ ≥ 1 || ndw * θ^(max_iter - iter) > κtol * (1 - θ)
        break
      end
    end

    @. w1 = w1 + dw1
    @. w2 = w2 + dw2
    @. w3 = w3 + dw3

    # transform `w` to `z`
    @. z1 = T11 * w1 + T12 * w2 + T13 * w3
    @. z2 = T21 * w1 + T22 * w2 + T23 * w3
    @. z3 = T31 * w1 +       w2           # T32 = 1, T33 = 0

    # check stopping criterion
    if iter == 1
      η = max(cache.ηold, eps(eltype(reltol)))^(0.8)
      do_newton = iszero(integrator.success_iter) || iter < min_iter || η * ndw > κtol
    else
      η = θ / (1 - θ) # calculated for possible early stopping
      do_newton = iter < min_iter || η * ndw > κtol
    end
  end
  cache.ηold = η
  integrator.force_stepfail = do_newton
  cache.nl_iters = iter
  do_newton && return

  @. u = uprev + z3

  if adaptive
    utilde = w2
    e1dt, e2dt, e3dt = e1/dt, e2/dt, e3/dt
    @. tmp = e1dt*z1 + e2dt*z2 + e3dt*z3
    mass_matrix != I && (mul!(w1, mass_matrix, tmp); copyto!(tmp, w1))
    @. utilde = integrator.fsalfirst + tmp
    alg.smooth_est && (linsolve1(vec(utilde), W1, vec(utilde), false); integrator.destats.nsolve += 1)
    # RadauIIA5 needs a transformed rtol and atol see
    # https://github.com/luchr/ODEInterface.jl/blob/0bd134a5a358c4bc13e0fb6a90e27e4ee79e0115/src/radau5.f#L399-L421
    calculate_residuals!(atmp, utilde, uprev, u, atol, rtol, internalnorm,t)
    integrator.EEst = internalnorm(atmp,t)

    if !(integrator.EEst < oneunit(integrator.EEst)) && integrator.iter == 1 || integrator.u_modified
      @. utilde = uprev + utilde
      f(fsallast, utilde, p, t)
      integrator.destats.nf += 1
      @. utilde = fsallast + tmp
      alg.smooth_est && (linsolve1(vec(utilde), W1, vec(utilde), false); integrator.destats.nsolve += 1)
      calculate_residuals!(atmp, utilde, uprev, u, atol, rtol, internalnorm,t)
      integrator.EEst = internalnorm(atmp,t)
    end
  end

  if integrator.EEst <= oneunit(integrator.EEst)
    cache.dtprev = dt
    if alg.extrapolant != :constant
      @. cache.cont1  = (z2 - z3)/c2m1
      @. tmp          = (z1 - z2)/c1mc2
      @. cache.cont2  = (tmp - cache.cont1) / c1m1
      @. cache.cont3  = cache.cont2 - (tmp - z1/c1)/c2
    end
  end

  f(fsallast, u, p, t+dt)
  integrator.destats.nf += 1
  return
end
