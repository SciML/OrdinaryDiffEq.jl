function do_newJ(integrator, alg, cache, repeat_step)::Bool # for FIRK
  integrator.iter <= 1 && return true
  repeat_step && return false
  first(islinearfunction(integrator)) && return false
  integrator.opts.adaptive || return true
  alg_can_repeat_jac(alg) || return true
  integrator.u_modified && return true
  # below is Newton specific logic, so we return non-Newton algs here
  alg isa NewtonAlgorithm || return true
  nlstatus = cache.status
  Int8(nlstatus) < 0 && return true
  # no reuse when the cutoff is 0
  fast_convergence_cutoff = alg.fast_convergence_cutoff
  iszero(fast_convergence_cutoff) && return true
  # reuse J when there is fast convergence
  fastconvergence = nlstatus === DiffEqBase.FastConvergence
  return !fastconvergence
end

function do_newW(integrator, nlsolver, new_jac, W_dt)::Bool # for FIRK
  nlsolver === nothing && return true
  new_jac && return true
  # reuse W when the change in stepsize is small enough
  dt = integrator.dt
  smallstepchange = abs((dt-W_dt)/W_dt) <= get_new_W_γdt_cutoff(nlsolver)
  return !smallstepchange
end

function initialize!(integrator,cache::RadauIIA3ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  nothing
end

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

function initialize!(integrator, cache::RadauIIA3Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
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
      @.. cache.rtol = reltol^(2/3) / 10
      @.. cache.atol = cache.rtol * (abstol / reltol)
    end
  end
  nothing
end

@muladd function perform_step!(integrator,cache::RadauIIA3ConstantCache)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack T11, T12, T21, T22, TI11, TI12, TI21, TI22 = cache.tab
  @unpack c1, c2, α, β, e1, e2 = cache.tab
  @unpack κ, cont1, cont2 = cache
  @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
  alg = unwrap_alg(integrator, true)
  @unpack maxiters = alg
  mass_matrix = integrator.f.mass_matrix

  # precalculations
  rtol = @. reltol^(2/3) / 10
  atol = @. rtol * (abstol / reltol)
  αdt, βdt = α/dt, β/dt
  J = calc_J(integrator,  cache)

  cache.dtprev = one(cache.dtprev)
  z1 = w1 = map(zero, u)
  z2 = w2 = map(zero, u)
  cache.cont1 = map(zero, u)
  cache.cont2 = map(zero, u)

  if u isa Number
    LU1 = -(αdt + βdt*im)*mass_matrix + J
  else
    LU1 = lu(-(αdt + βdt*im)*mass_matrix + J)
  end
  integrator.destats.nw += 1

  # Newton iteration
  local ndw, ff1, ff2
  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  fail_convergence = true
  iter = 0
  while iter < maxiters
    iter += 1
    integrator.destats.nnonliniter += 1
    # evaluate function
    ff1 = f(uprev+z1, p, t+c1*dt)
    ff2 = f(uprev+z2, p, t+c2*dt)
    integrator.destats.nf += 2

    fw1 = @. TI11 * ff1 + TI12 * ff2
    fw2 = @. TI21 * ff1 + TI22 * ff2

    if mass_matrix isa UniformScaling
      Mw1 = @. mass_matrix.λ * w1
      Mw2 = @. mass_matrix.λ * w2
    else
      Mw1 = mass_matrix*w1
      Mw2 = mass_matrix*w2
    end

    rhs1 = @. fw1 - αdt*Mw1 + βdt*Mw2
    rhs2 = @. fw2 - βdt*Mw1 - αdt*Mw2
    dw12 = LU1 \ (@. rhs1 + rhs2*im)
    integrator.destats.nsolve += 1
    dw1 = real(dw12)
    dw2 = imag(dw12)

    # compute norm of residuals
    iter > 1 && (ndwprev = ndw)
    atmp1 = calculate_residuals(dw1, uprev, u, atol, rtol, internalnorm, t)
    atmp2 = calculate_residuals(dw2, uprev, u, atol, rtol, internalnorm, t)
    ndw = internalnorm(atmp1, t) + internalnorm(atmp2, t)

    # check divergence (not in initial step)
    if iter > 1
      θ = ndw / ndwprev
      ( diverge = θ > 1 ) && ( cache.status = DiffEqBase.Divergence )
      ( veryslowconvergence = ndw * θ^(maxiters - iter) > κ * (1 - θ) ) && ( cache.status = DiffEqBase.VerySlowConvergence )
      if diverge || veryslowconvergence
        break
      end
    end

    w1 = @. w1 - dw1
    w2 = @. w2 - dw2

    # transform `w` to `z`
    z1 = @. T11 * w1 + T12 * w2
    z2 = @. T21 * w1 + T22 * w2

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
      # Newton method converges
      cache.status = η < alg.fast_convergence_cutoff ? DiffEqBase.FastConvergence : DiffEqBase.Convergence
      fail_convergence = false
      break
    end
  end

  cache.ηold = η
  cache.iter = iter

  u = @. uprev + z2

  if adaptive
    utilde = @. dt*(e1*ff1 + e2*ff2)
    atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm, t)
    integrator.EEst = internalnorm(atmp, t)
  end

  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
  return
end

@muladd function perform_step!(integrator, cache::RadauIIA3Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,fsallast,fsalfirst = integrator
  @unpack T11, T12, T21, T22, TI11, TI12, TI21, TI22 = cache.tab
  @unpack c1, c2, α, β, e1, e2 = cache.tab
  @unpack κ, cont1, cont2 = cache
  @unpack z1, z2, w1, w2,
          dw12, cubuff,
          k, k2, fw1, fw2,
          J, W1,
          tmp, atmp, jac_config, rtol, atol = cache
  @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
  alg = unwrap_alg(integrator, true)
  @unpack maxiters = alg
  mass_matrix = integrator.f.mass_matrix
  # precalculations
  αdt, βdt = α/dt, β/dt
  (new_jac = do_newJ(integrator, alg, cache, repeat_step)) && (calc_J!(J, integrator, cache); cache.W_γdt = dt)
  if (new_W = do_newW(integrator, alg, new_jac, cache.W_γdt))
      @inbounds for II in CartesianIndices(J)
        W1[II] = -(αdt + βdt*im) * mass_matrix[Tuple(II)...] + J[II]
      end
      integrator.destats.nw += 1
  end

  #better initial guess
  uzero = zero(eltype(z1))
  @. z1 = uzero
  @. z2 = uzero
  @. w1 = uzero
  @. w2 = uzero
  @. cache.cont1 = uzero
  @. cache.cont2 = uzero

  # Newton iteration
  local ndw
  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  fail_convergence = true
  iter = 0
  while iter < maxiters
    iter += 1
    integrator.destats.nnonliniter += 1
    # evaluate function
    @. tmp = uprev + z1
    f(fsallast, tmp, p, t+c1*dt)
    @. tmp = uprev + z2
    f(k2, tmp, p, t+c2*dt)
    integrator.destats.nf += 2

    @. fw1 = TI11 * fsallast + TI12 * k2
    @. fw2 = TI21 * fsallast + TI22 * k2

    if mass_matrix === I
      Mw1 = w1
      Mw2 = w2
    elseif mass_matrix isa UniformScaling
      mul!(z1, mass_matrix.λ, w1)
      mul!(z2, mass_matrix.λ, w2)
      Mw1 = z1
      Mw2 = z2
    else
      mul!(z1, mass_matrix, w1)
      mul!(z2, mass_matrix, w2)
      Mw1 = z1
      Mw2 = z2
    end

    @. cubuff = complex(fw1 - αdt*Mw1 + βdt*Mw2, fw2 - βdt*Mw1 - αdt*Mw2)
    needfactor = iter==1

    linsolve = cache.linsolve
    if needfactor
      linsolve = LinearSolve.set_A(linsolve,W1)
    end
    linres = dolinsolve(integrator, linsolve; b = _vec(cubuff), u = _vec(dw12))
    cache.linsolve = linres.cache

    integrator.destats.nsolve += 1
    dw1 = real(dw12)
    dw2 = imag(dw12)

    # compute norm of residuals
    iter > 1 && (ndwprev = ndw)
    calculate_residuals!(atmp, dw1, uprev, u, atol, rtol, internalnorm, t)
    ndw1 = internalnorm(atmp, t)
    calculate_residuals!(atmp, dw2, uprev, u, atol, rtol, internalnorm, t)
    ndw2 = internalnorm(atmp, t)
    ndw = ndw1 + ndw2

    # check divergence (not in initial step)
    if iter > 1
      θ = ndw / ndwprev
      ( diverge = θ > 2 ) && ( cache.status = DiffEqBase.Divergence )
      if diverge
        break
      end
    end

    @. w1 = w1 - dw1
    @. w2 = w2 - dw2

    # transform `w` to `z`
    @. z1 = T11 * w1 + T12 * w2
    @. z2 = T21 * w1 + T22 * w2

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
      # Newton method converges
      cache.status = η < alg.fast_convergence_cutoff ? DiffEqBase.FastConvergence : DiffEqBase.Convergence
      fail_convergence = false
      break
    end
  end
  if fail_convergence
    integrator.force_stepfail = true
    integrator.destats.nnonlinconvfail += 1
    return
  end
  cache.ηold = η
  cache.iter = iter

  @. u = uprev + z2
  if adaptive
    utilde = w2
    @. utilde = dt*(e1*fsallast + e2*k2)
    calculate_residuals!(atmp, utilde, uprev, u, atol, rtol, internalnorm, t)
    integrator.EEst = internalnorm(atmp, t)
  end

  f(fsallast, u, p, t+dt)
  integrator.destats.nf += 1
  return
end

@muladd function perform_step!(integrator, cache::RadauIIA5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack T11, T12, T13, T21, T22, T23, T31, TI11, TI12, TI13, TI21, TI22, TI23, TI31, TI32, TI33 = cache.tab
  @unpack c1, c2, γ, α, β, e1, e2, e3 = cache.tab
  @unpack κ, cont1, cont2, cont3 = cache
  @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
  alg = unwrap_alg(integrator, true)
  @unpack maxiters = alg
  mass_matrix = integrator.f.mass_matrix

  # precalculations
  rtol = @.. reltol^(2/3) / 10
  atol = @.. rtol * (abstol / reltol)
  c1m1 = c1-1
  c2m1 = c2-1
  c1mc2= c1-c2
  γdt, αdt, βdt = γ/dt, α/dt, β/dt
  J = calc_J(integrator,  cache)
  if u isa Number
    LU1 = -γdt*mass_matrix + J
    LU2 = -(αdt + βdt*im)*mass_matrix + J
  else
    LU1 = lu(-γdt*mass_matrix + J)
    LU2 = lu(-(αdt + βdt*im)*mass_matrix + J)
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
    z1  = @.. c1′ * (cont1 + (c1′-c2m1) * (cont2 + (c1′-c1m1) * cont3))
    z2  = @.. c2′ * (cont1 + (c2′-c2m1) * (cont2 + (c2′-c1m1) * cont3))
    z3  = @.. c3′ * (cont1 + (c3′-c2m1) * (cont2 + (c3′-c1m1) * cont3))
    w1  = @.. TI11 * z1 + TI12 * z2 + TI13 * z3
    w2  = @.. TI21 * z1 + TI22 * z2 + TI23 * z3
    w3  = @.. TI31 * z1 + TI32 * z2 + TI33 * z3
  end

  # Newton iteration
  local ndw
  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  fail_convergence = true
  iter = 0
  while iter < maxiters
    iter += 1
    integrator.destats.nnonliniter += 1

    # evaluate function
    ff1 = f(uprev+z1, p, t+c1*dt)
    ff2 = f(uprev+z2, p, t+c2*dt)
    ff3 = f(uprev+z3, p, t+   dt) # c3 = 1
    integrator.destats.nf += 3

    fw1 = @.. TI11 * ff1 + TI12 * ff2 + TI13 * ff3
    fw2 = @.. TI21 * ff1 + TI22 * ff2 + TI23 * ff3
    fw3 = @.. TI31 * ff1 + TI32 * ff2 + TI33 * ff3

    if mass_matrix isa UniformScaling # `UniformScaling` doesn't play nicely with broadcast
      Mw1 = @.. mass_matrix.λ * w1
      Mw2 = @.. mass_matrix.λ * w2
      Mw3 = @.. mass_matrix.λ * w3
    else
      Mw1 = mass_matrix*w1
      Mw2 = mass_matrix*w2
      Mw3 = mass_matrix*w3
    end

    rhs1 = @.. fw1 - γdt*Mw1
    rhs2 = @.. fw2 - αdt*Mw2 + βdt*Mw3
    rhs3 = @.. fw3 - βdt*Mw2 - αdt*Mw3
    dw1 = LU1 \ rhs1
    dw23 = LU2 \ (@.. rhs2 + rhs3*im)
    integrator.destats.nsolve += 2
    dw2 = real(dw23)
    dw3 = imag(dw23)

    # compute norm of residuals
    iter > 1 && (ndwprev = ndw)
    atmp1 = calculate_residuals(dw1, uprev, u, atol, rtol, internalnorm, t)
    atmp2 = calculate_residuals(dw2, uprev, u, atol, rtol, internalnorm, t)
    atmp3 = calculate_residuals(dw3, uprev, u, atol, rtol, internalnorm, t)
    ndw = internalnorm(atmp1, t) + internalnorm(atmp2, t) + internalnorm(atmp3, t)

    # check divergence (not in initial step)
    if iter > 1
      θ = ndw / ndwprev
      ( diverge = θ > 1 ) && ( cache.status = DiffEqBase.Divergence )
      ( veryslowconvergence = ndw * θ^(maxiters - iter) > κ * (1 - θ) ) && ( cache.status = DiffEqBase.VerySlowConvergence )
      if diverge || veryslowconvergence
        break
      end
    end

    w1 = @.. w1 - dw1
    w2 = @.. w2 - dw2
    w3 = @.. w3 - dw3

    # transform `w` to `z`
    z1 = @.. T11 * w1 + T12 * w2 + T13 * w3
    z2 = @.. T21 * w1 + T22 * w2 + T23 * w3
    z3 = @.. T31 * w1 +       w2           # T32 = 1, T33 = 0

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
      # Newton method converges
      cache.status = η < alg.fast_convergence_cutoff ? DiffEqBase.FastConvergence : DiffEqBase.Convergence
      fail_convergence = false
      break
    end
  end
  if fail_convergence
    integrator.force_stepfail = true
    integrator.destats.nnonlinconvfail += 1
    return
  end
  cache.ηold = η
  cache.iter = iter

  u = @.. uprev + z3

  if adaptive
    e1dt, e2dt, e3dt = e1/dt, e2/dt, e3/dt
    tmp = @.. e1dt*z1 + e2dt*z2 + e3dt*z3
    mass_matrix != I && (tmp = mass_matrix*tmp)
    utilde = @.. integrator.fsalfirst + tmp
    alg.smooth_est && (utilde = LU1 \ utilde; integrator.destats.nsolve += 1)
    # RadauIIA5 needs a transformed rtol and atol see
    # https://github.com/luchr/ODEInterface.jl/blob/0bd134a5a358c4bc13e0fb6a90e27e4ee79e0115/src/radau5.f#L399-L421
    atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm, t)
    integrator.EEst = internalnorm(atmp, t)

    if !(integrator.EEst < oneunit(integrator.EEst)) && integrator.iter == 1 || integrator.u_modified
      f0 = f(uprev .+ utilde, p, t)
      integrator.destats.nf += 1
      utilde = @.. f0 + tmp
      alg.smooth_est && (utilde = LU1 \ utilde; integrator.destats.nsolve += 1)
      atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm, t)
      integrator.EEst = internalnorm(atmp, t)
    end
  end

  if integrator.EEst <= oneunit(integrator.EEst)
    cache.dtprev = dt
    if alg.extrapolant != :constant
      cache.cont1  = @.. (z2 - z3)/c2m1
      tmp          = @.. (z1 - z2)/c1mc2
      cache.cont2  = @.. (tmp - cache.cont1) / c1m1
      cache.cont3  = @.. cache.cont2 - (tmp - z1/c1)/c2
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
  @unpack κ, cont1, cont2, cont3 = cache
  @unpack z1, z2, z3, w1, w2, w3,
          dw1, ubuff, dw23, cubuff,
          k, k2, k3, fw1, fw2, fw3,
          J, W1, W2,
          tmp, atmp, jac_config, linsolve1, linsolve2, rtol, atol = cache
  @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
  alg = unwrap_alg(integrator, true)
  @unpack maxiters = alg
  mass_matrix = integrator.f.mass_matrix

  # precalculations
  c1m1 = c1-1
  c2m1 = c2-1
  c1mc2= c1-c2
  γdt, αdt, βdt = γ/dt, α/dt, β/dt
  (new_jac = do_newJ(integrator, alg, cache, repeat_step)) && (calc_J!(J, integrator, cache); cache.W_γdt = dt)
  if (new_W = do_newW(integrator, alg, new_jac, cache.W_γdt))
    @inbounds for II in CartesianIndices(J)
      W1[II] = -γdt * mass_matrix[Tuple(II)...] + J[II]
      W2[II] = -(αdt + βdt*im) * mass_matrix[Tuple(II)...] + J[II]
    end
    integrator.destats.nw += 1
  end

  # TODO better initial guess
  if integrator.iter == 1 || integrator.u_modified || alg.extrapolant == :constant
    cache.dtprev = one(cache.dtprev)
    uzero = zero(eltype(u))
    @.. z1 = uzero
    @.. z2 = uzero
    @.. z3 = uzero
    @.. w1 = uzero
    @.. w2 = uzero
    @.. w3 = uzero
    @.. cache.cont1 = uzero
    @.. cache.cont2 = uzero
    @.. cache.cont3 = uzero
  else
    c3′ = dt/cache.dtprev
    c1′ = c1*c3′
    c2′ = c2*c3′
    @.. z1  = c1′ * (cont1 + (c1′-c2m1) * (cont2 + (c1′-c1m1) * cont3))
    @.. z2  = c2′ * (cont1 + (c2′-c2m1) * (cont2 + (c2′-c1m1) * cont3))
    @.. z3  = c3′ * (cont1 + (c3′-c2m1) * (cont2 + (c3′-c1m1) * cont3))
    @.. w1  = TI11 * z1 + TI12 * z2 + TI13 * z3
    @.. w2  = TI21 * z1 + TI22 * z2 + TI23 * z3
    @.. w3  = TI31 * z1 + TI32 * z2 + TI33 * z3
  end

  # Newton iteration
  local ndw
  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  fail_convergence = true
  iter = 0
  while iter < maxiters
    iter += 1
    integrator.destats.nnonliniter += 1

    # evaluate function
    @.. tmp = uprev + z1
    f(fsallast, tmp, p, t+c1*dt)
    @.. tmp = uprev + z2
    f(k2, tmp, p, t+c2*dt)
    @.. tmp = uprev + z3
    f(k3, tmp, p, t+   dt) # c3 = 1
    integrator.destats.nf += 3

    @.. fw1 = TI11 * fsallast + TI12 * k2 + TI13 * k3
    @.. fw2 = TI21 * fsallast + TI22 * k2 + TI23 * k3
    @.. fw3 = TI31 * fsallast + TI32 * k2 + TI33 * k3

    if mass_matrix === I
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

    @.. ubuff = fw1 - γdt*Mw1
    needfactor = iter==1 && new_W

    linsolve1 = cache.linsolve1
    if needfactor
      linsolve1 = LinearSolve.set_A(linsolve1,W1)
    end
    linres1 = dolinsolve(integrator, linsolve1; b = _vec(ubuff), u = _vec(dw1))
    cache.linsolve1 = linres1.cache

    @.. cubuff = complex(fw2 - αdt*Mw2 + βdt*Mw3, fw3 - βdt*Mw2 - αdt*Mw3)

    linsolve2 = cache.linsolve2
    if needfactor
      linsolve2 = LinearSolve.set_A(linsolve2,W2)
    end
    linres2 = dolinsolve(integrator, linsolve2; b = _vec(cubuff), u = _vec(dw23))
    cache.linsolve2 = linres2.cache

    integrator.destats.nsolve += 2
    dw2 = z2; dw3 = z3
    @.. dw2 = real(dw23)
    @.. dw3 = imag(dw23)

    # compute norm of residuals
    iter > 1 && (ndwprev = ndw)
    calculate_residuals!(atmp, dw1, uprev, u, atol, rtol, internalnorm, t)
    ndw1 = internalnorm(atmp, t)
    calculate_residuals!(atmp, dw2, uprev, u, atol, rtol, internalnorm, t)
    ndw2 = internalnorm(atmp, t)
    calculate_residuals!(atmp, dw3, uprev, u, atol, rtol, internalnorm, t)
    ndw3 = internalnorm(atmp, t)
    ndw = ndw1 + ndw2 + ndw3

    # check divergence (not in initial step)
    if iter > 1
      θ = ndw / ndwprev
      ( diverge = θ > 1 ) && ( cache.status = DiffEqBase.Divergence )
      ( veryslowconvergence = ndw * θ^(maxiters - iter) > κ * (1 - θ) ) && ( cache.status = DiffEqBase.VerySlowConvergence )
      if diverge || veryslowconvergence
        break
      end
    end

    @.. w1 = w1 - dw1
    @.. w2 = w2 - dw2
    @.. w3 = w3 - dw3

    # transform `w` to `z`
    @.. z1 = T11 * w1 + T12 * w2 + T13 * w3
    @.. z2 = T21 * w1 + T22 * w2 + T23 * w3
    @.. z3 = T31 * w1 +       w2           # T32 = 1, T33 = 0

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
      # Newton method converges
      cache.status = η < alg.fast_convergence_cutoff ? DiffEqBase.FastConvergence : DiffEqBase.Convergence
      fail_convergence = false
      break
    end
  end
  if fail_convergence
    integrator.force_stepfail = true
    integrator.destats.nnonlinconvfail += 1
    return
  end
  cache.ηold = η
  cache.iter = iter

  @.. u = uprev + z3

  if adaptive
    utilde = w2
    e1dt, e2dt, e3dt = e1/dt, e2/dt, e3/dt
    @.. tmp = e1dt*z1 + e2dt*z2 + e3dt*z3
    mass_matrix != I && (mul!(w1, mass_matrix, tmp); copyto!(tmp, w1))
    @.. ubuff = integrator.fsalfirst + tmp

    if alg.smooth_est
      linres1 = dolinsolve(integrator, linsolve1; b = _vec(ubuff), u = _vec(utilde))
      cache.linsolve1 = linres1.cache
      integrator.destats.nsolve += 1
    end

    # RadauIIA5 needs a transformed rtol and atol see
    # https://github.com/luchr/ODEInterface.jl/blob/0bd134a5a358c4bc13e0fb6a90e27e4ee79e0115/src/radau5.f#L399-L421
    calculate_residuals!(atmp, utilde, uprev, u, atol, rtol, internalnorm, t)
    integrator.EEst = internalnorm(atmp, t)

    if !(integrator.EEst < oneunit(integrator.EEst)) && integrator.iter == 1 || integrator.u_modified
      @.. utilde = uprev + utilde
      f(fsallast, utilde, p, t)
      integrator.destats.nf += 1
      @.. ubuff = fsallast + tmp

      if alg.smooth_est
        linres1 = dolinsolve(integrator, linsolve1; b = _vec(ubuff), u = _vec(utilde))
        cache.linsolve1 = linres1.cache
        integrator.destats.nsolve += 1
      end

      calculate_residuals!(atmp, utilde, uprev, u, atol, rtol, internalnorm, t)
      integrator.EEst = internalnorm(atmp, t)
    end
  end

  if integrator.EEst <= oneunit(integrator.EEst)
    cache.dtprev = dt
    if alg.extrapolant != :constant
      @.. cache.cont1  = (z2 - z3)/c2m1
      @.. tmp          = (z1 - z2)/c1mc2
      @.. cache.cont2  = (tmp - cache.cont1) / c1m1
      @.. cache.cont3  = cache.cont2 - (tmp - z1/c1)/c2
    end
  end

  f(fsallast, u, p, t+dt)
  integrator.destats.nf += 1
  return
end
