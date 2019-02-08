@muladd function (S::NLAnderson{false})(integrator)
  nlcache = S.cache
  @unpack t,dt,uprev,u,p = integrator
  @unpack z,tmp,κ,tol,c,γ,max_iter,min_iter,zs,gs,residuals = nlcache
  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)
  # precalculations
  κtol = κ*tol

  # initial step of fixed point iteration with Anderson acceleration
  iter = 1
  tstep = t + c*dt
  u = @. tmp + γ*z
  if mass_matrix == I
    z₊ = dt .* f(u, p, tstep)
  else
    mz = _reshape(mass_matrix * _vec(z), axes(z))
    z₊ = dt .* f(u, p, tstep) .- mz .+ z
  end
  zs[1] = z
  gs[1] = z₊

  # circularly shift zs and gs and update zs[1]
  zptr = zs[S.n+1]
  gptr = gs[S.n+1]
  for t in S.n:-1:1
    zs[t + 1] = zs[t]
    gs[t + 1] = gs[t]
  end
  zs[1] = z₊

  # compute initial values for early stopping criterion
  ndz = integrator.opts.internalnorm(z₊ .- z,t)

  # update solution
  z = z₊

  # check stopping criterion for initial step
  η = nlcache.ηold
  do_anderson = true # TODO: this makes `min_iter` ≥ 2

  # fixed point iteration with Anderson acceleration
  while do_anderson && iter < max_iter
    # compute next iterate
    iter += 1
    u = @. tmp + γ*z
    if mass_matrix == I
      z₊ = dt .* f(u, p, tstep)
    else
      mz = _reshape(mass_matrix * _vec(z), axes(z))
      z₊ = dt .* f(u, p, tstep) .- mz .+ z
    end
    gs[1] = z₊

    mk = min(S.n, iter-1)
    ztmp = zs[1] - gs[1]
    for i in 2:(mk + 1)
      residuals[:, i - 1] .= _vec(gs[i]) .- _vec(zs[i]) .+ _vec(ztmp)
    end
    resqr = qr!(view(residuals, :, 1:mk), Val(true))
    alphas = ztmp isa Number ? resqr \ fill(ztmp, 1, 1) : resqr \ _vec(ztmp)
    for i in 1:mk
        z₊ = @. z₊ + alphas[i] * (gs[i+1] - gs[1])
    end

    # circularly shift zs and gs and update zs[1]
    zptr = zs[S.n+1]
    gptr = gs[S.n+1]
    for t in S.n:-1:1
      zs[t + 1] = zs[t]
      gs[t + 1] = gs[t]
    end
    zs[1] = z₊

    # check early stopping criterion
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(z₊ .- z,t)
    θ = ndz/ndzprev
    if θ ≥ 1 || ndz * θ^(max_iter - iter) > κtol * (1 - θ)
      break
    end

    # update solution
    z = z₊

    # check stopping criterion
    η = θ / (1 - θ) # calculated for possible early stopping
    do_anderson = iter < min_iter || η * ndz > κtol
  end

  integrator.force_stepfail = do_anderson
  z, η, iter, do_anderson
end

@muladd function (S::NLAnderson{true})(integrator)
  nlcache = S.cache
  @unpack t,dt,uprev,u,p = integrator
  @unpack z,z₊,dz,ztmp,tmp,κ,tol,k,c,γ,min_iter,max_iter,zs,gs,alphas,residuals = nlcache
  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)
  # precalculations
  κtol = κ*tol
  vecztmp = vec(ztmp); vecz = vec(z); vecz₊ = vec(z₊)

  # initial step of fixed point iteration with Anderson acceleration
  iter = 1
  tstep = t + c*dt
  @. u = tmp + γ*z
  f(k, u, p, tstep)
  if mass_matrix == I
    @. z₊ = dt*k
  else
    @. z₊ = dt*k + z
    mul!(vecztmp, mass_matrix, vecz)
    @. z₊ -= ztmp
  end
  zs[1] .= vecz
  gs[1] .= vecz₊

  # circularly shift zs and gs and update zs[1]
  zptr = zs[S.n+1]
  gptr = gs[S.n+1]
  for t in S.n:-1:1
    zs[t + 1] = zs[t]
    gs[t + 1] = gs[t]
  end
  zs[1] = zptr # required to not update zs[2] implicitly in the next line
  zs[1] .= vecz₊
  gs[1] = gptr

  # compute initial values for early stopping criterion
  @. dz = z₊ - z
  ndz = integrator.opts.internalnorm(dz,t)

  # update solution
  z .= z₊

  # check stopping criterion for initial step
  η = nlcache.ηold
  do_anderson = true # TODO: this makes `min_iter` ≥ 2

  # fixed point iteration with Anderson acceleration
  while do_anderson && iter < max_iter
    # compute next iterate
    iter += 1
    @. u = tmp + γ*z
    f(k, u, p, tstep)
    if mass_matrix == I
      @. z₊ = dt*k
    else
      @. z₊ = dt*k + z
      mul!(vecztmp, mass_matrix, vecz)
      @. z₊ -= ztmp
    end
    gs[1] .= vecz₊

    mk = min(S.n, iter-1)
    @. vecztmp = zs[1] - gs[1]
    for i in 1:mk
      @. residuals[:, i] = gs[i + 1] - zs[i + 1] + vecztmp
    end
    resqr = qr!(view(residuals, :, 1:mk), Val(true))
    ldiv!(view(alphas, 1:mk), resqr, vecztmp)
    for i in 1:mk
        @. vecz₊ += alphas[i] * (gs[i + 1] - gs[1])
    end

    # circularly shift zs and gs and update zs[1]
    zptr = zs[S.n+1]
    gptr = gs[S.n+1]
    for t in S.n:-1:1
      zs[t + 1] = zs[t]
      gs[t + 1] = gs[t]
    end
    zs[1] = zptr # required to not update zs[2] implicitly in the next line
    zs[1] .= vecz₊
    gs[1] = gptr

    # check early stopping criterion
    ndzprev = ndz
    @. dz = z₊ - z
    ndz = integrator.opts.internalnorm(dz,t)
    θ = ndz/ndzprev
    if θ ≥ 1 || ndz * θ^(max_iter - iter) > κtol * (1 - θ)
      break
    end

    # update solution
    @. z = z₊

    # check stopping criterion
    η = θ / (1 - θ) # calculated for possible early stopping
    do_anderson = iter < min_iter || η * ndz > κtol
  end

  integrator.force_stepfail = do_anderson
  z, η, iter, do_anderson
end
