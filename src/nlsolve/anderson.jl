function (S::NLAnderson{false,<:NLSolverCache})(integrator)
  nlcache = S.cache
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,tmp,W,κ,tol,c,γ,max_iter,min_iter = nlcache
  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
  else
    f = integrator.f
  end
  # precalculations
  κtol = κ*tol

  zs = zeros(S.n+1)
  gs = zeros(S.n+1)
  residuals = zeros(S.n+1)
  alphas = zeros(S.n)
  # initial step of NLAnderson iteration
  zs[1] = z
  iter = 1
  tstep = t + c*dt
  u = tmp + γ*z
  z₊ = dt*f(u, p, tstep)
  gs[1] = z₊
  dz = z₊ - z
  ndz = integrator.opts.internalnorm(dz)
  for t = (S.n+1):-1:2
    zs[t] = zs[t-1]
    gs[t] = gs[t-1]
  end
  # zs = circshift(zs, 1)
  # gs = circshift(gs, 1)
  z = z₊
  zs[1] = z
  η = nlcache.ηold
  do_anderson = true

  # anderson acceleration for fixed point iteration
  fail_convergence = false
  while (do_anderson || iter < min_iter) && iter < max_iter
    iter += 1
    u = tmp + γ*z
    z₊ = dt*f(u, p, tstep)
    gs[1] = z₊
    
    mk = min(S.n, iter-1)
    residuals[1:mk] = (gs[2:mk+1] .- zs[2:mk+1]) .- (gs[1] - zs[1])
    alphas[1:mk] .= residuals[1:mk] \ [(zs[1] - gs[1])]
    for i = 1:mk
        z₊ += alphas[i]*(gs[i+1] - gs[1])
    end
    for t = (S.n+1):-1:2
      zs[t] = zs[t-1]
      gs[t] = gs[t-1]
    end
    # zs = circshift(zs, 1)
    # gs = circshift(gs, 1)
    zs[1] = z₊
    ndzprev = ndz
    dz = z₊ - z
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(max_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_anderson = (η*ndz > κtol)
    z = z₊
  end

  if (iter >= max_iter && do_anderson) || fail_convergence
    integrator.force_stepfail = true
    return (z, η, iter, true)
  end

  return (z, η, iter, false)
end

function (S::NLAnderson{true})(integrator)
  nlcache = S.cache
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,z₊,b,dz,tmp,κ,tol,k,c,γ,min_iter,max_iter = nlcache
  ztmp = b
  mass_matrix = integrator.f.mass_matrix
  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
  else
    f = integrator.f
  end
  # precalculations
  κtol = κ*tol

  zs = fill(Any[], S.n+1)
  gs = fill(Any[], S.n+1)
  residuals = zeros(length(vec(z)), S.n+1)
  for e = 1:1:(S.n+1)
    zs[e] = zeros(length(vec(z)))
    gs[e] = zeros(length(vec(z)))
  end
  alphas = zeros(S.n)
  # initial step of NLAnderson iteration
  zs[1] .= vec(z)
  iter = 1
  tstep = t + c*dt
  @. u = tmp + γ*z
  # z₊ = dt*f(u, p, tstep)
  f(k, u, p, tstep)
  if mass_matrix == I
    @. z₊ = dt*k
  else
    @. z₊ = dt*k + z
    mul!(ztmp, mass_matrix, z)
    @. z₊ -= ztmp
  end
  gs[1] .= vec(z₊)
  @. dz = z₊ - z
  ndz = integrator.opts.internalnorm(dz)
  for t = (S.n+1):-1:2
    zs[t] = zs[t-1]
    gs[t] = gs[t-1]
  end
  # zs = circshift(zs, 1)
  # gs = circshift(gs, 1)
  z .= z₊
  zs[1] .= vec(z)
  η = nlcache.ηold
  do_anderson = true
  # anderson acceleration for fixed point iteration
  fail_convergence = false
  while (do_anderson || iter < min_iter) && iter < max_iter
    iter += 1
    @. u = tmp + γ*z
    f(k, u, p, tstep)
    if mass_matrix == I
      @. z₊ = dt*k
    else
      @. z₊ = dt*k + z
      mul!(ztmp, mass_matrix, z)
      @. z₊ -= ztmp
    end
    gs[1] .= vec(z₊)

    mk = min(S.n, iter-1)
    tmp = (gs[1] - zs[1])
    for i in 2:mk+1
      residuals[:,i-1] .= (gs[i] .- zs[i]) .- (gs[1] - zs[1])
    end
    alphas[1:mk] .= residuals[:,1:mk] \ (zs[1] .- gs[1])
    vecz₊ = vec(z₊)
    for i = 1:mk
        vecz₊ .+= alphas[i].*(gs[i+1] .- gs[1])
    end
    for t = (S.n+1):-1:2
      zs[t] = zs[t-1]
      gs[t] = gs[t-1]
    end
    # zs = circshift(zs, 1)
    # gs = circshift(gs, 1)
    zs[1] .= vec(z₊)
    ndzprev = ndz
    @. dz = z₊ - z
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(max_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_anderson = (η*ndz > κtol)
    @. z = z₊
  end

  if (iter >= max_iter && do_anderson) || fail_convergence
    integrator.force_stepfail = true
    return (z, η, iter, true)
  end

  return (z, η, iter, false)
end