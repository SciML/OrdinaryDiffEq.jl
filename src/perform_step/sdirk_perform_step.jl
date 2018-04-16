function initialize!(integrator, cache::ImplicitEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ImplicitEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,κ,tol = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  # precalculations
  κtol = κ*tol

  # calculate W
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*J
  end

  # initial guess
  if alg.extrapolant == :linear
    z = dt*integrator.fsalfirst
  else # :constant
    z = zero(u)
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  u = uprev + z
  b = dt*f(u, p, tstep) - z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = uprev + z
    b = dt*f(u, p, tstep) - z
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z = z + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = uprev + z

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive && integrator.success_iter > 0
    # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
    # use 2nd divided differences (DD) a la SPICE and Shampine

    # TODO: check numerical stability
    uprev2 = integrator.uprev2
    tprev = integrator.tprev

    dt1 = dt*(t+dt-tprev)
    dt2 = (t-tprev)*(t+dt-tprev)
    c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
    r = c*dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

    tmp = r*abs((u - uprev)/dt1 - (uprev - uprev2)/dt2)
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  else
    integrator.EEst = 1
  end

  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end#

function initialize!(integrator, cache::ImplicitEulerCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::ImplicitEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,du1,dz,z,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  mass_matrix = integrator.sol.prob.mass_matrix
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  # precalculations
  κtol = κ*tol

  new_W = calc_W!(integrator, cache, dt, repeat_step)

  # initial guess
  if alg.extrapolant == :linear
    @. z = dt*integrator.fsalfirst
  else # :constant
    z .= zero(u)
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. u = uprev + z
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
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = uprev + z
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
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = uprev + z

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive && integrator.success_iter > 0
    # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
    # use 2nd divided differences (DD) a la SPICE and Shampine

    # TODO: check numerical stability
    uprev2 = integrator.uprev2
    tprev = integrator.tprev

    dt1 = dt*(t+dt-tprev)
    dt2 = (t-tprev)*(t+dt-tprev)
    c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
    r = c*dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

    @. tmp = r*abs((u - uprev)/dt1 - (uprev - uprev2)/dt2)
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  else
    integrator.EEst = 1
  end

  f(integrator.fsallast,u,p,t+dt)
end

function initialize!(integrator, cache::ImplicitMidpointConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ImplicitMidpointConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,κ,tol = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

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
  if alg.extrapolant == :linear
    z = dt*integrator.fsalfirst
  else # :constant
    z = zero(u)
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dto2
  # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
  u = uprev + z/2
  b = dt*f(u, p, tstep) - z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
    u = uprev + z/2
    b = dt*f(u, p, tstep) - z
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z = z + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = uprev + z

  cache.ηold = η
  cache.newton_iters = iter

  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::ImplicitMidpointCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::ImplicitMidpointCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,du1,dz,z,k,b,J,W,jac_config,tmp,κ,tol = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  mass_matrix = integrator.sol.prob.mass_matrix

  # precalculations
  κtol = κ*tol

  dto2 = dt/2

  new_W = calc_W!(integrator, cache, dto2, repeat_step)

  # initial guess
  if alg.extrapolant == :linear
    @. z = dt*integrator.fsalfirst
  else # :constant
    z .= zero(u)
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dto2
  # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
  @. u = uprev + z/2
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
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
    @. u = uprev + z/2
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
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = uprev + z

  cache.ηold = η
  cache.newton_iters = iter

  f(integrator.fsallast,u,p,t+dt)
end

function initialize!(integrator, cache::TrapezoidConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::TrapezoidConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,κ,tol = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

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
  zprev = dt*integrator.fsalfirst
  z = zprev # Constant extrapolation

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = uprev + dto2*integrator.fsalfirst
  u = tmp + z/2
  b = dt*f(u, p, tstep) - z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + z/2
    b = dt*f(u, p, tstep) - z
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z = z + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = tmp + z/2

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    if integrator.iter > 2
      # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
      # use 3rd divided differences (DD) a la SPICE and Shampine

      # TODO: check numerical stability
      uprev2 = integrator.uprev2
      tprev = integrator.tprev
      uprev3 = cache.uprev3
      tprev2 = cache.tprev2

      dt1 = dt*(t+dt-tprev)
      dt2 = (t-tprev)*(t+dt-tprev)
      dt3 = (t-tprev)*(t-tprev2)
      dt4 = (tprev-tprev2)*(t-tprev2)
      dt5 = t+dt-tprev2
      c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
      r = c*dt^3/2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

      # tmp = r*abs(((u - uprev)/dt1 - (uprev - uprev2)/dt2) - ((uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4)/dt5)
      DD31 = (u - uprev)/dt1 - (uprev - uprev2)/dt2
      DD30 = (uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4
      tmp = r*abs((DD31 - DD30)/dt5)
      atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp)
      if integrator.EEst <= 1
        cache.uprev3 = uprev2
        cache.tprev2 = tprev
      end
    elseif integrator.success_iter > 0
      integrator.EEst = 1
      cache.uprev3 = integrator.uprev2
      cache.tprev2 = integrator.tprev
    else
      integrator.EEst = 1
    end
  end

  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::TrapezoidCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::TrapezoidCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,du1,dz,z,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  mass_matrix = integrator.sol.prob.mass_matrix

  # precalculations
  κtol = κ*tol

  dto2 = dt/2

  new_W = calc_W!(integrator, cache, dto2, repeat_step)

  # initial guess
  @. z = dt*integrator.fsalfirst

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. tmp = uprev + dto2*integrator.fsalfirst
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
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
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
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + z/2

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    if integrator.iter > 2
      # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
      # use 3rd divided differences (DD) a la SPICE and Shampine

      # TODO: check numerical stability
      uprev2 = integrator.uprev2
      tprev = integrator.tprev
      uprev3 = cache.uprev3
      tprev2 = cache.tprev2

      dt1 = dt*(t+dt-tprev)
      dt2 = (t-tprev)*(t+dt-tprev)
      dt3 = (t-tprev)*(t-tprev2)
      dt4 = (tprev-tprev2)*(t-tprev2)
      dt5 = t+dt-tprev2
      c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
      r = c*dt^3/2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

      # @. tmp = r*abs(((u - uprev)/dt1 - (uprev - uprev2)/dt2) - ((uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4)/dt5)
      @inbounds for i in eachindex(u)
        DD31 = (u[i] - uprev[i])/dt1 - (uprev[i] - uprev2[i])/dt2
        DD30 = (uprev[i] - uprev2[i])/dt3 - (uprev2[i] - uprev3[i])/dt4
        tmp[i] = r*abs((DD31 - DD30)/dt5)
      end
      calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp)
      if integrator.EEst <= 1
        copy!(cache.uprev3,uprev2)
        cache.tprev2 = tprev
      end
    elseif integrator.success_iter > 0
      integrator.EEst = 1
      copy!(cache.uprev3,integrator.uprev2)
      cache.tprev2 = integrator.tprev
    else
      integrator.EEst = 1
    end
  end

  f(integrator.fsallast,u,p,t+dt)
end

function initialize!(integrator, cache::TRBDF2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::TRBDF2ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,d,ω,btilde1,btilde2,btilde3,α1,α2 = cache.tab
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  # precalculations
  κtol = κ*tol

  ddt = d*dt

  # calculate W
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - ddt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - ddt*J
  end

  # FSAL
  zprev = dt*integrator.fsalfirst

  ##### Solve Trapezoid Step

  # TODO: Add extrapolation
  zᵧ = zprev

  # initial step of Newton iteration
  iter = 1
  tstep = t + γ*dt
  tmp = uprev + d*zprev
  u = tmp + d*zᵧ
  b = dt*f(u, p, tstep) - zᵧ
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  zᵧ = zᵧ + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + d*zᵧ
    b = dt*f(u, p, tstep) - zᵧ
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    zᵧ = zᵧ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve BDF2 Step

  ### Initial Guess From Shampine
  z = α1*zprev + α2*zᵧ

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = uprev + ω*zprev + ω*zᵧ
  u = tmp + d*z
  b = dt*f(u, p, tstep) - z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + d*z
    b = dt*f(u, p, tstep) - z
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z = z + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = tmp + d*z

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = btilde1*zprev + btilde2*zᵧ + btilde3*z
    if alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = z./dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::TRBDF2Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::TRBDF2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,du1,dz,zprev,zᵧ,z,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,d,ω,btilde1,btilde2,btilde3,α1,α2 = cache.tab
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  # precalculations
  κtol = κ*tol

  ddt = d*dt

  new_W = calc_W!(integrator, cache, ddt, repeat_step)

  # FSAL
  @. zprev = dt*integrator.fsalfirst

  ##### Solve Trapezoid Step

  # TODO: Add extrapolation
  @. zᵧ = zprev

  # initial step of Newton iteration
  iter = 1
  tstep = t + γ*dt
  @. tmp = uprev + d*zprev
  @. u = tmp + d*zᵧ
  f(k, u, p, tstep)
  @. b = dt*k - zᵧ
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  zᵧ .+= dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + d*zᵧ
    f(k, u, p, tstep)
    @. b = dt*k - zᵧ
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    zᵧ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve BDF2 Step

  ### Initial Guess From Shampine
  @. z = α1*zprev + α2*zᵧ

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. tmp = uprev + ω*zprev + ω*zᵧ
  @. u = tmp + d*z
  f(k, u, p, tstep)
  @. b = dt*k - z
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + d*z
    f(k, u, p, tstep)
    @. b = dt*k - z
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + d*z

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    @. dz = btilde1*zprev + btilde2*zᵧ + btilde3*z
    if alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        cache.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z/dt
end

function initialize!(integrator, cache::SDIRK2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::SDIRK2ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,κ,tol = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  # precalculations
  κtol = κ*tol

  # calculate W
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*J
  end

  # initial guess
  if integrator.success_iter > 0 && !integrator.reeval_fsal && alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
    z₁ = u - uprev
  elseif alg.extrapolant == :linear
    z₁ = dt*integrator.fsalfirst
  else
    z₁ = zero(u)
  end

  ##### Step 1

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  u = uprev + z₁
  b = dt*f(u, p, tstep) - z₁
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₁ = z₁ + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = uprev + z₁
    b = dt*f(u, p, tstep) - z₁
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₁ = z₁ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 2

  ### Initial Guess Is α₁ = c₂/γ, c₂ = 0 => z₂ = α₁z₁ = 0
  z₂ = zero(u)

  # initial step of Newton iteration
  iter = 1
  tmp = uprev - z₁
  u = tmp .+ z₂
  b = dt*f(u, p, t) - z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp .+ z₂
    b = dt*f(u, p, t) - z₂
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ = z₂ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = uprev + z₁/2 + z₂/2

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = z₁/2 - z₂/2
    if alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = f(u, p, t)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::SDIRK2Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::SDIRK2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,du1,dz,z₁,z₂,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  # precalculations
  κtol = κ*tol

  new_W = calc_W!(integrator, cache, dt, repeat_step)

  # initial guess
  if integrator.success_iter > 0 && !integrator.reeval_fsal && alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
    @. z₁ = u - uprev
  elseif alg.extrapolant == :linear
    @. z₁ = dt*integrator.fsalfirst
  else
    z₁ .= zero(u)
  end

  ##### Step 1

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. u = uprev + z₁
  f(k, u, p, tstep)
  @. b = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₁ .+= dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = uprev + z₁
    f(k, u, p, tstep)
    @. b = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₁ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 2

  ### Initial Guess Is α₁ = c₂/γ, c₂ = 0 => z₂ = α₁z₁ = 0
  z₂ .= zero(u)

  # initial step of Newton iteration
  iter = 1
  @. tmp = uprev - z₁
  @. u = tmp + z₂
  f(k, u, p, t)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    @. u = tmp + z₂
    f(k, u, p, t)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = uprev + z₁/2 + z₂/2

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    @. dz = z₁/2 - z₂/2
    if alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        cache.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  f(integrator.fsallast,u,p,t)
end

function initialize!(integrator, cache::SSPSDIRK2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::SSPSDIRK2ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,κ,tol = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  γ = eltype(u)(1//4)
  c2 = typeof(t)(3//4)

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  # calculate W
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - γdt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - γdt*J
  end

  # initial guess
  if integrator.success_iter > 0 && !integrator.reeval_fsal && alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
    z₁ = u - uprev
  elseif alg.extrapolant == :linear
    z₁ = dt*integrator.fsalfirst
  else
    z₁ = zero(u)
  end

  ##### Step 1

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  u = uprev + γ*z₁
  b = dt*f(u, p, tstep) - z₁
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₁ = z₁ + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = uprev + γ*z₁
    b = dt*f(u, p, tstep) - z₁
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₁ = z₁ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 2

  ### Initial Guess Is α₁ = c₂/γ
  z₂ = c2/γ

  # initial step of Newton iteration
  iter = 1
  tmp = uprev + z₁/2
  u = tmp + z₂/4
  b = dt*f(u, p, t) - z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + z₂/4
    b = dt*f(u, p, t) - z₂
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ = z₂ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = tmp + z₂/2

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  integrator.fsallast = f(u, p, t)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::SSPSDIRK2Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::SSPSDIRK2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,du1,dz,z₁,z₂,k,b,J,W,jac_config,tmp,κ,tol = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  γ = eltype(u)(1//4)
  c2 = typeof(t)(3//4)

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  new_W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  if integrator.success_iter > 0 && !integrator.reeval_fsal && alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
    @. z₁ = u - uprev
  elseif alg.extrapolant == :linear
    @. z₁ = dt*integrator.fsalfirst
  else
    z₁ .= zero(u)
  end

  ##### Step 1

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. u = uprev + γ*z₁
  f(k, u, p, tstep)
  @. b = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₁ .+= dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁
    f(k, u, p, tstep)
    @. b = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₁ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 2

  ### Initial Guess Is α₁ = c₂/γ
  @. z₂ = c2/γ

  iter = 1
  @. tmp = uprev + z₁/2
  @. u = tmp + z₂/4
  f(k, u, p, t)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    @. u = tmp + z₂/4
    f(k, u, p, t)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + z₂/2

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  f(integrator.fsallast,u,p,t)
end

function initialize!(integrator, cache::Cash4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Cash4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack b1hat1,b2hat1,b3hat1,b4hat1,b1hat2,b2hat2,b3hat2,b4hat2 = cache.tab
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  # calculate W
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - γdt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - γdt*J
  end

  ##### Step 1

  # TODO: Add extrapolation for guess
  z₁ = zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + γdt
  u = uprev + γ*z₁
  b = dt*f(u, p, tstep) - z₁
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₁ = z₁ + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = uprev + γ*z₁
    b = dt*f(u, p, tstep) - z₁
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₁ = z₁ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ##### Step 2

  # TODO: Add extrapolation for guess
  z₂ = zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + c2*dt
  tmp = uprev + a21*z₁
  u = tmp .+ γ*z₂
  b = dt*f(u, p, tstep) - z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₂
    b = dt*f(u, p, tstep) - z₂
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ = z₂ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # Guess starts from z₁
  z₃ = z₁

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  tmp = uprev + a31*z₁ + a32*z₂
  u = tmp + γ*z₃
  b = dt*f(u, p, tstep) - z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₃
    b = dt*f(u, p, tstep) - z₃
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₃ = z₃ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  z₄ = z₃

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃
  u = tmp + γ*z₄
  b = dt*f(u, p, tstep) - z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₄
    b = dt*f(u, p, tstep) - z₄
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₄ = z₄ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat2 for prediction
  z₅ = b1hat2*z₁ + b2hat2*z₂ + b3hat2*z₃ + b4hat2*z₄

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  u = tmp + γ*z₅
  b = dt*f(u, p, tstep) - z₅
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₅ = z₅ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₅
    b = dt*f(u, p, tstep) - z₅
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₅ = z₅ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = tmp + γ*z₅

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    if alg.embedding == 3
      btilde1 = b1hat2-a51; btilde2 = b2hat2-a52;
      btilde3 = b3hat2-a53; btilde4 = b4hat2-a54; btilde5 = -γ
    else
      btilde1 = b1hat1-a51; btilde2 = b2hat1-a52;
      btilde3 = b3hat1-a53; btilde4 = b4hat1-a54; btilde5 = -γ
    end

    tmp = btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅
    if alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = z₅./dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::Cash4Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::Cash4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,z₅,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack b1hat1,b2hat1,b3hat1,b4hat1,b1hat2,b2hat2,b3hat2,b4hat2 = cache.tab
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  new_W = calc_W!(integrator, cache, γdt, repeat_step)

  ##### Step 1

  # TODO: Add extrapolation for guess
  z₁ .= zero(z₁)

  # initial step of Newton iteration
  iter = 1
  tstep = t + γdt
  @. u = uprev + γ*z₁
  f(k, u, p, tstep)
  @. b = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₁ .+= dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁
    f(k, u, p, tstep)
    @. b = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₁ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ##### Step 2

  # TODO: Add extrapolation for guess
  z₂ .= zero(z₂)

  # initial step of Newton iteration
  iter = 1
  tstep = t + c2*dt
  @. tmp = uprev + a21*z₁
  @. u = tmp + γ*z₂
  f(k, u, p, tstep)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₂
    f(k, u, p, tstep)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # Guess starts from z₁
  @. z₃ = z₁

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  @. tmp = uprev + a31*z₁ + a32*z₂
  @. u = tmp + γ*z₃
  f(k, u, p, tstep)
  @. b = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₃ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₃
    f(k, u, p, tstep)
    @. b = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₃ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = z₃

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  @. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃
  @. u = tmp + γ*z₄
  f(k, u, p, tstep)
  @. b = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₄ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₄
    f(k, u, p, tstep)
    @. b = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₄ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use constant z prediction
  @. z₅ = b1hat2*z₁ + b2hat2*z₂ + b3hat2*z₃ + b4hat2*z₄

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  @. u = tmp + γ*z₅
  f(k, u, p, tstep)
  @. b = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₅ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₅
    f(k, u, p, tstep)
    @. b = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₅ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + γ*z₅

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    if alg.embedding == 3
      btilde1 = b1hat2-a51; btilde2 = b2hat2-a52;
      btilde3 = b3hat2-a53; btilde4 = b4hat2-a54; btilde5 = -γ
    else
      btilde1 = b1hat1-a51; btilde2 = b2hat1-a52;
      btilde3 = b3hat1-a53; btilde4 = b4hat1-a54; btilde5 = -γ
    end

    @. dz = btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅
    if alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        cache.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z₅/dt
end

function initialize!(integrator, cache::Hairer4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Hairer4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack α21,α31,α32,α41,α43 = cache.tab
  @unpack bhat1,bhat2,bhat3,bhat4,btilde1,btilde2,btilde3,btilde4,btilde5 = cache.tab
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  # calculate W
  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - γdt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - γdt*J
  end

  # TODO: Add extrapolation for guess
  z₁ = zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + γdt
  u = uprev + γ*z₁
  b = dt*f(u, p, tstep) - z₁
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₁ = z₁ + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = uprev + γ*z₁
    b = dt*f(u, p, tstep) - z₁
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₁ = z₁ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ##### Step 2

  z₂ = α21*z₁

  # initial step of Newton iteration
  iter = 1
  tstep = t + c2*dt
  tmp = uprev + a21*z₁
  u = tmp + γ*z₂
  b = dt*f(u, p, tstep) - z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₂
    b = dt*f(u, p, tstep) - z₂
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ = z₂ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  z₃ = α31*z₁ + α32*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  tmp = uprev + a31*z₁ + a32*z₂
  u = tmp + γ*z₃
  b = dt*f(u, p, tstep) - z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₃
    b = dt*f(u, p, tstep) - z₃
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₃ = z₃ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  z₄ = α41*z₁ + α43*z₃

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃
  u = tmp + γ*z₄
  b = dt*f(u, p, tstep) - z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₄
    b = dt*f(u, p, tstep) - z₄
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₄ = z₄ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat2 for prediction
  z₅ = bhat1*z₁ + bhat2*z₂ + bhat3*z₃ + bhat4*z₄

  # initial step of Newton iteration
  iter = 1
  tstep = t+dt
  tmp = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  u = tmp + γ*z₅
  b = dt*f(u, p, tstep) - z₅
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₅ = z₅ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₅
    b = dt*f(u, p, tstep) - z₅
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₅ = z₅ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = tmp + γ*z₅

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅
    if alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = z₅./dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::Hairer4Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::Hairer4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,z₅,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack α21,α31,α32,α41,α43 = cache.tab
  @unpack bhat1,bhat2,bhat3,bhat4,btilde1,btilde2,btilde3,btilde4,btilde5 = cache.tab
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  new_W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  if integrator.success_iter > 0 && !integrator.reeval_fsal && alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
    @. z₁ = u - uprev
  elseif alg.extrapolant == :linear
    @. z₁ = dt*integrator.fsalfirst
  else
    z₁ .= zero(0)
  end

  ##### Step 1

  # initial step of Newton iteration
  iter = 1
  tstep = t + γdt
  @. u = uprev + γ*z₁
  f(k, u, p, tstep)
  @. b = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₁ .+= dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁
    f(k, u, p, tstep)
    @. b = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₁ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ##### Step 2

  @. z₂ = α21*z₁

  # initial step of Newton iteration
  iter = 1
  tstep = t + c2*dt
  @. tmp = uprev + a21*z₁
  @. u = tmp + γ*z₂
  f(k, u, p, tstep)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₂
    f(k, u, p, tstep)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  @. z₃ = α31*z₁ + α32*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  @. tmp = uprev + a31*z₁ + a32*z₂
  @. u = tmp + γ*z₃
  f(k, u, p, tstep)
  @. b = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₃ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₃
    f(k, u, p, tstep)
    @. b = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₃ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = α41*z₁ + α43*z₃

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  @. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃
  @. u = tmp + γ*z₄
  f(k, u, p, tstep)
  @. b = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₄ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₄
    f(k, u, p, tstep)
    @. b = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₄ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat prediction
  @. z₅ = bhat1*z₁ + bhat2*z₂ + bhat3*z₃ + bhat4*z₄

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. tmp = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  @. u = tmp + γ*z₅
  f(k, u, p, tstep)
  @. b = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₅ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₅
    f(k, u, p, tstep)
    @. b = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₅ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + γ*z₅

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    # @. dz = btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅
    @tight_loop_macros for i in eachindex(u)
      dz[i] = btilde1*z₁[i] + btilde2*z₂[i] + btilde3*z₃[i] + btilde4*z₄[i] + btilde5*z₅[i]
    end
    if alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        cache.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z₅/dt
end
