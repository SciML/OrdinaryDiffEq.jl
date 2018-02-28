function initialize!(integrator, cache::CNABConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::CNABConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,κ,tol = cache
  @unpack a0,a1 = cache

  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
    f2 = integrator.f.f2
  else
    f = integrator.f
  end

  # precalculations
  κtol = κ*tol

  dto2 = dt/2

  # calculate W
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dto2*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dto2*J
  end

  # initial guess
  zprev = 0
  if typeof(integrator.f) <: SplitFunction
    zprev = dt.*f(uprev, p, t)
  else
    zprev = dt*integrator.fsalfirst
  end
  z = zprev # Constant extrapolation

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = uprev + dto2*integrator.fsalfirst

  if typeof(integrator.f) <: SplitFunction
    # This assumes the implicit part is cheaper than the explicit part
    y0 = dt*integrator.fsalfirst - z
    tmp += a0*y0
  end

  u = tmp + z/2
  b = dt*f(u, p, tstep) - z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = tmp + z/2
    b = dt*f(u, p, tstep) - z
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

  if typeof(integrator.f) <: SplitFunction
    u = tmp + z/2
    y1 = dt*f2(u,p,tstep)
    tmp += a1*y1 + z/2
    u = tmp
  else
    u = tmp + z/2
  end




  cache.ηold = η
  cache.newton_iters = iter
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::CNABCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::CNABCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,du1,dz,z,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
    f2 = integrator.f.f2
  else
    f = integrator.f
  end

  # precalculations
  κtol = κ*tol

  dto2 = dt/2

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},W,uprev,p,dto2,t) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac}, J, uprev, p, t)
      else
        uf.t = t
        jacobian!(J, uf, uprev, du1, integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(typeof(integrator.t)))
      new_W = true
      for j in 1:length(u), i in 1:length(u)
        @inbounds W[i,j] = mass_matrix[i,j]-dto2*J[i,j]
      end
    else
      new_W = false
    end
  end

  if typeof(integrator.f) <: SplitFunction
    # Explicit tableau is not FSAL
    # Make this not compute on repeat
    if !repeat_step && !integrator.last_stepfail
      f(z, integrator.uprev, p, integrator.t)
      z .*= dt
    end
  else
    # FSAL Step 1
    @. z = dt*integrator.fsalfirst
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. tmp = uprev + dto2*integrator.fsalfirst

  if typeof(integrator.f) <: SplitFunction
    # This assumes the implicit part is cheaper than the explicit part
    @. y0 = dt*integrator.fsalfirst - z
    @. tmp += a0*y0
  end

  @. u = tmp + z/2
  f(k, u, p, tstep)
  if mass_matrix == I
    @. b = dt*k - z
  else
    A_mul_B!(vec(b),mass_matrix,vec(z))
    @. b = dt*k - b
  end
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + z/2
    f(k, u, p, tstep)
    if mass_matrix == I
      @. b = dt*k - z
    else
      A_mul_B!(vec(b),mass_matrix,vec(z))
      @. b = dt*k - b
    end
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

  

  cache.ηold = η
  cache.newton_iters = iter

  if typeof(integrator.f) <: SplitFunction
    @. u = tmp + z/2
    f2(y2, u, p, tstep); y2 .*= dt
      for i in eachindex(tmp)
        @inbounds tmp[i] += uprev[i] + a1*y2[i]
      end
      @. u = tmp
  else
    @. u = tmp + z/2
  end

  f(integrator.fsallast,u,p,t+dt)
end