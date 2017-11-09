function initialize!(integrator, cache::ImplicitEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ImplicitEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache

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
  if integrator.alg.extrapolant == :linear
    z = dt.*integrator.fsalfirst
  else # :constant
    z = zero(u)
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  u = uprev .+ z
  b = dt.*f(tstep,u) .- z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = uprev .+ z
    b = dt.*f(tstep,u) .- z
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
    z = z .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = uprev .+ z

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

    tmp = @. r*abs((u - uprev)/dt1 - (uprev - uprev2)/dt2)
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  else
    integrator.EEst = 1
  end

  integrator.fsallast = f(t+dt,u)
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
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::ImplicitEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  # precalculations
  κtol = κ*tol

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,dt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      for j in 1:length(u), i in 1:length(u)
        @inbounds W[i,j] = mass_matrix[i,j]-dt*J[i,j]
      end
    else
      new_W = false
    end
  end

  # initial guess
  if integrator.alg.extrapolant == :linear
    @. z = dt*integrator.fsalfirst
  else # :constant
    z .= zero(u)
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. u = uprev + z
  f(tstep,u,k)
  if mass_matrix == I
    @. b = dt*k - z
  else
    A_mul_B!(vec(b),mass_matrix,vec(z))
    @. b = dt*k - b
  end
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + z
    f(tstep,u,k)
    if mass_matrix == I
      @. b = dt*k - z
    else
      A_mul_B!(vec(b),mass_matrix,vec(z))
      @. b = dt*k - b
    end
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  else
    integrator.EEst = 1
  end

  f(t+dt,u,integrator.fsallast)
end

function initialize!(integrator, cache::ImplicitMidpointConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ImplicitMidpointConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache

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
  if integrator.alg.extrapolant == :linear
    z = dt.*integrator.fsalfirst
  else # :constant
    z = zero(u)
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dto2
  # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
  u = @. uprev + z/2
  b = dt.*f(tstep,u) .- z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
    u = @. uprev + z/2
    b = dt.*f(tstep,u) .- z
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
    z = z .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. uprev + z

  cache.ηold = η
  cache.newton_iters = iter

  integrator.fsallast = f(t+dt,u)
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
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::ImplicitMidpointCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z,k,b,J,W,jac_config,tmp,κ,tol = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  # precalculations
  κtol = κ*tol

  dto2 = dt/2

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,dto2,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-dto2*J[i,j]
      end
    else
      new_W = false
    end
  end

  # initial guess
  if integrator.alg.extrapolant == :linear
    @. z = dt*integrator.fsalfirst
  else # :constant
    z .= zero(u)
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dto2
  # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
  @. u = uprev + z/2
  f(tstep,u,k)
  if mass_matrix == I
    @. b = dt*k - z
  else
    A_mul_B!(vec(b),mass_matrix,vec(z))
    @. b = dt*k - b
  end
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
    @. u = uprev + z/2
    f(tstep,u,k)
    if mass_matrix == I
      @. b = dt*k - z
    else
      A_mul_B!(vec(b),mass_matrix,vec(z))
      @. b = dt*k - b
    end
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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

  @. u = uprev + z

  cache.ηold = η
  cache.newton_iters = iter

  f(t+dt,u,integrator.fsallast)
end

function initialize!(integrator, cache::TrapezoidConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::TrapezoidConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache

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
  zprev = dt.*integrator.fsalfirst
  z = zprev # Constant extrapolation

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = @. uprev + dto2*integrator.fsalfirst
  u = @. tmp + z/2
  b = dt.*f(tstep,u) .- z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + z/2
    b = dt.*f(tstep,u) .- z
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
    z = z .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. tmp + z/2

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

      # tmp = @. r*abs(((u - uprev)/dt1 - (uprev - uprev2)/dt2) - ((uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4)/dt5)
      DD31 = @. (u - uprev)/dt1 - (uprev - uprev2)/dt2
      DD30 = @. (uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4
      tmp = @. r*abs((DD31 - DD30)/dt5)
      atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
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

  integrator.fsallast = f(t+dt,u)
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
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::TrapezoidCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  # precalculations
  κtol = κ*tol

  dto2 = dt/2

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,dto2,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        if alg_autodiff(integrator.alg)
          ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
        else
          Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
        end
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      for j in 1:length(u), i in 1:length(u)
        @inbounds W[i,j] = mass_matrix[i,j]-dto2*J[i,j]
      end
    else
      new_W = false
    end
  end

  # initial guess
  @. z = dt*integrator.fsalfirst

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. tmp = uprev + dto2*integrator.fsalfirst
  @. u = tmp + z/2
  f(tstep,u,k)
  if mass_matrix == I
    @. b = dt*k - z
  else
    A_mul_B!(vec(b),mass_matrix,vec(z))
    @. b = dt*k - b
  end
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + z/2
    f(tstep,u,k)
    if mass_matrix == I
      @. b = dt*k - z
    else
      A_mul_B!(vec(b),mass_matrix,vec(z))
      @. b = dt*k - b
    end
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
      calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
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

  f(t+dt,u,integrator.fsallast)
end

function initialize!(integrator, cache::TRBDF2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::TRBDF2ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,d,ω,btilde1,btilde2,btilde3,α1,α2 = cache.tab

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
  zprev = dt.*integrator.fsalfirst

  ##### Solve Trapezoid Step

  # TODO: Add extrapolation
  zᵧ = zprev

  # initial step of Newton iteration
  iter = 1
  tstep = t + γ*dt
  tmp = @. uprev + d*zprev
  u = @. tmp + d*zᵧ
  b = dt.*f(tstep,u) .- zᵧ
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  zᵧ = zᵧ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + d*zᵧ
    b = dt.*f(tstep,u) .- zᵧ
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
    zᵧ = zᵧ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve BDF2 Step

  ### Initial Guess From Shampine
  z = @. α1*zprev + α2*zᵧ

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = @. uprev + ω*zprev + ω*zᵧ
  u = @. tmp + d*z
  b = dt.*f(tstep,u) .- z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + d*z
    b = dt.*f(tstep,u) .- z
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
    z = z .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. tmp + d*z

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. btilde1*zprev + btilde2*zᵧ + btilde3*z
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
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
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::TRBDF2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,zprev,zᵧ,z,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,d,ω,btilde1,btilde2,btilde3,α1,α2 = cache.tab

  # precalculations
  κtol = κ*tol

  ddt = d*dt

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,ddt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      mass_matrix = integrator.sol.prob.mass_matrix
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-ddt*J[i,j]
      end
    else
      new_W = false
    end
  end

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
  f(tstep,u,k)
  @. b = dt*k - zᵧ
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  zᵧ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + d*zᵧ
    f(tstep,u,k)
    @. b = dt*k - zᵧ
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    zᵧ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve BDF2 Step

  ### Initial Guess From Shampine
  z = @. α1*zprev + α2*zᵧ

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. tmp = uprev + ω*zprev + ω*zᵧ
  @. u = tmp + d*z
  f(tstep,u,k)
  @. b = dt*k - z
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + d*z
    f(tstep,u,k)
    @. b = dt*k - z
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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

  @. u = tmp + d*z

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    @. dz = btilde1*zprev + btilde2*zᵧ + btilde3*z
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        integrator.alg.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z/dt
end

function initialize!(integrator, cache::SDIRK2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::SDIRK2ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache

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
  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
    z₁ = u .- uprev
  elseif integrator.alg.extrapolant == :linear
    z₁ = dt.*integrator.fsalfirst
  else
    z₁ = zero(u)
  end

  ##### Step 1

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  u = @. uprev + z₁
  b = dt.*f(tstep,u) .- z₁
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₁ = z₁ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + z₁
    b = dt.*f(tstep,u) .- z₁
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
    z₁ = z₁ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 2

  ### Initial Guess Is α₁ = c₂/γ, c₂ = 0 => z₂ = α₁z₁ = 0
  z₂ = zero(u)

  # initial step of Newton iteration
  iter = 1
  tmp = @. uprev - z₁
  u = tmp .+ z₂
  b = dt.*f(t,u) .- z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = tmp .+ z₂
    b = dt.*f(t,u) .- z₂
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
    z₂ = z₂ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. uprev + z₁/2 + z₂/2

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. z₁/2 - z₂/2
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = f(t,u)
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
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::SDIRK2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z₁,z₂,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache

  # precalculations
  κtol = κ*tol

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,dt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      mass_matrix = integrator.sol.prob.mass_matrix
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-dt*J[i,j]
      end
    else
      new_W = false
    end
  end

  # initial guess
  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
    @. z₁ = u - uprev
  elseif integrator.alg.extrapolant == :linear
    @. z₁ = dt*integrator.fsalfirst
  else
    z₁ .= zero(u)
  end

  ##### Step 1

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. u = uprev + z₁
  f(tstep,u,k)
  @. b = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₁ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + z₁
    f(tstep,u,k)
    @. b = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₁ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(t,u,k)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    @. u = tmp + z₂
    f(t,u,k)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₂ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = uprev + z₁/2 + z₂/2

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    @. dz = z₁/2 - z₂/2
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        integrator.alg.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  f(t,u,integrator.fsallast)
end

function initialize!(integrator, cache::SSPSDIRK2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::SSPSDIRK2ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache

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
  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
    z₁ = @. u - uprev
  elseif integrator.alg.extrapolant == :linear
    z₁ = @. dt*integrator.fsalfirst
  else
    z₁ = zero(u)
  end

  ##### Step 1

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  u = @. uprev + γ*z₁
  b = dt.*f(tstep,u) .- z₁
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₁ = z₁ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + γ*z₁
    b = dt.*f(tstep,u) .- z₁
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
    z₁ = z₁ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 2

  ### Initial Guess Is α₁ = c₂/γ
  z₂ = c2/γ

  # initial step of Newton iteration
  iter = 1
  tmp = @. uprev + z₁/2
  u = @. tmp + z₂/4
  b = dt.*f(t,u) .- z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + z₂/4
    b = dt.*f(t,u) .- z₂
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
    z₂ = z₂ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. tmp + z₂/2

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  integrator.fsallast = f(t,u)
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
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::SSPSDIRK2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z₁,z₂,k,b,J,W,jac_config,tmp,κ,tol = cache

  γ = eltype(u)(1//4)
  c2 = typeof(t)(3//4)

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,γdt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      mass_matrix = integrator.sol.prob.mass_matrix
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
      end
    else
      new_W = false
    end
  end

  # initial guess
  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
    @. z₁ = u - uprev
  elseif integrator.alg.extrapolant == :linear
    @. z₁ = dt*integrator.fsalfirst
  else
    z₁ .= zero(u)
  end

  ##### Step 1

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. u = uprev + γ*z₁
  f(tstep,u,k)
  @. b = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₁ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁
    f(tstep,u,k)
    @. b = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₁ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 2

  ### Initial Guess Is α₁ = c₂/γ
  @. z₂ = c2/γ

  iter = 1
  @. tmp = uprev + z₁/2
  @. u = tmp + z₂/4
  f(t,u,k)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    @. u = tmp + z₂/4
    f(t,u,k)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₂ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + z₂/2

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  f(t,u,integrator.fsallast)
end

function initialize!(integrator, cache::Union{Kvaerno3ConstantCache,KenCarp3ConstantCache})
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Union{Kvaerno3ConstantCache,KenCarp3ConstantCache}, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,α32 = cache.tab

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

  # FSAL Step 1
  z₁ = dt.*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation for guess
  z₂ = z₁

  iter = 1
  tstep = t + 2*γdt
  tmp = @. uprev + γ*z₁
  u = @. tmp + γ*z₂
  b = dt.*f(tstep,u) .- z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₂
    b = dt.*f(tstep,u) .- z₂
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
    z₂ = z₂ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # Guess is from Hermite derivative on z₁ and z₂
  z₃ = @. α31*z₁ + α32*z₂

  iter = 1
  tstep = t + c3*dt
  tmp = @. uprev + a31*z₁ + a32*z₂
  u = @. tmp + γ*z₃
  b = dt.*f(tstep,u) .- z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₃
    b = dt.*f(tstep,u) .- z₃
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
    z₃ = z₃ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  if typeof(cache) <: Kvaerno3ConstantCache
    z₄ = @. a31*z₁ + a32*z₂ + γ*z₃ # use yhat as prediction
  elseif typeof(cache) <: KenCarp3ConstantCache
    @unpack α41,α42 = cache.tab
    z₄ = @. α41*z₁ + α42*z₂
  end

  iter = 1
  tstep = t + dt
  tmp = @. uprev + a41*z₁ + a42*z₂ + a43*z₃
  u = @. tmp + γ*z₄
  b = dt.*f(tstep,u) .- z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₄
    b = dt.*f(tstep,u) .- z₄
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
    z₄ = z₄ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. tmp + γ*z₄

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = z₄./dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::Union{Kvaerno3Cache,KenCarp3Cache})
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::Union{Kvaerno3Cache,KenCarp3Cache}, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,α32 = cache.tab

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,γdt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      mass_matrix = integrator.sol.prob.mass_matrix
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
      end
    else
      new_W = false
    end
  end

  # FSAL Step 1
  @. z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation for guess
  @. z₂ = z₁

  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt
  @. tmp = uprev + γ*z₁
  @. u = tmp + γ*z₂
  f(tstep,u,k)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₂
    f(tstep,u,k)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₂ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # Guess is from Hermite derivative on z₁ and z₂
  @. z₃ = α31*z₁ + α32*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  @. tmp = uprev + a31*z₁ + a32*z₂
  @. u = tmp + γ*z₃
  f(tstep,u,k)
  @. b = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₃ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₃
    f(tstep,u,k)
    @. b = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₃ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  if typeof(cache) <: Kvaerno3Cache
    @. z₄ = a31*z₁ + a32*z₂ + γ*z₃ # use yhat as prediction
  elseif typeof(cache) <: KenCarp3Cache
    @unpack α41,α42 = cache.tab
    @. z₄ = α41*z₁ + α42*z₂
  end

  iter = 1
  tstep = t + dt
  @. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃
  @. u = tmp + γ*z₄
  f(tstep,u,k)
  @. b = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₄ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₄
    f(tstep,u,k)
    @. b = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₄ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + γ*z₄

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    @. dz = btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        integrator.alg.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z₄/dt
end

function initialize!(integrator, cache::Cash4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Cash4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack b1hat1,b2hat1,b3hat1,b4hat1,b1hat2,b2hat2,b3hat2,b4hat2 = cache.tab

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
  u = @. uprev + γ*z₁
  b = dt.*f(tstep,u) .- z₁
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₁ = z₁ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + γ*z₁
    b = dt.*f(tstep,u) .- z₁
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
    z₁ = z₁ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ##### Step 2

  # TODO: Add extrapolation for guess
  z₂ = zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + c2*dt
  tmp = @. uprev + a21*z₁
  u = @. tmp .+ γ*z₂
  b = dt.*f(tstep,u) .- z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₂
    b = dt.*f(tstep,u) .- z₂
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
    z₂ = z₂ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # Guess starts from z₁
  z₃ = z₁

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  tmp = @. uprev + a31*z₁ + a32*z₂
  u = @. tmp + γ*z₃
  b = dt.*f(tstep,u) .- z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₃
    b = dt.*f(tstep,u) .- z₃
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
    z₃ = z₃ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  z₄ = z₃

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  tmp = @. uprev + a41*z₁ + a42*z₂ + a43*z₃
  u = @. tmp + γ*z₄
  b = dt.*f(tstep,u) .- z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₄
    b = dt.*f(tstep,u) .- z₄
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
    z₄ = z₄ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat2 for prediction
  z₅ = @. b1hat2*z₁ + b2hat2*z₂ + b3hat2*z₃ + b4hat2*z₄

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  u = @. tmp + γ*z₅
  b = dt.*f(tstep,u) .- z₅
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₅ = z₅ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₅
    b = dt.*f(tstep,u) .- z₅
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
    z₅ = z₅ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. tmp + γ*z₅

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    if integrator.alg.embedding == 3
      btilde1 = b1hat2-a51; btilde2 = b2hat2-a52;
      btilde3 = b3hat2-a53; btilde4 = b4hat2-a54; btilde5 = -γ
    else
      btilde1 = b1hat1-a51; btilde2 = b2hat1-a52;
      btilde3 = b3hat1-a53; btilde4 = b4hat1-a54; btilde5 = -γ
    end

    tmp = @. btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
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
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::Cash4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,z₅,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack b1hat1,b2hat1,b3hat1,b4hat1,b1hat2,b2hat2,b3hat2,b4hat2 = cache.tab

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,γdt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      mass_matrix = integrator.sol.prob.mass_matrix
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##### Step 1

  # TODO: Add extrapolation for guess
  z₁ .= zero(z₁)

  # initial step of Newton iteration
  iter = 1
  tstep = t + γdt
  @. u = uprev + γ*z₁
  f(tstep,u,k)
  @. b = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₁ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁
    f(tstep,u,k)
    @. b = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₁ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(tstep,u,k)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₂
    f(tstep,u,k)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₂ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(tstep,u,k)
  @. b = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₃ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₃
    f(tstep,u,k)
    @. b = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₃ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(tstep,u,k)
  @. b = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₄ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₄
    f(tstep,u,k)
    @. b = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₄ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use constant z prediction
  @. z₅ = b1hat2*z₁ + b2hat2*z₂ + b3hat2*z₃ + b4hat2*z₄

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  @. u = tmp + γ*z₅
  f(tstep,u,k)
  @. b = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₅ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₅
    f(tstep,u,k)
    @. b = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₅ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + γ*z₅

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    if integrator.alg.embedding == 3
      btilde1 = b1hat2-a51; btilde2 = b2hat2-a52;
      btilde3 = b3hat2-a53; btilde4 = b4hat2-a54; btilde5 = -γ
    else
      btilde1 = b1hat1-a51; btilde2 = b2hat1-a52;
      btilde3 = b3hat1-a53; btilde4 = b4hat1-a54; btilde5 = -γ
    end

    @. dz = btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        integrator.alg.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z₅/dt
end

function initialize!(integrator, cache::Hairer4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Hairer4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack α21,α31,α32,α41,α43 = cache.tab
  @unpack bhat1,bhat2,bhat3,bhat4,btilde1,btilde2,btilde3,btilde4,btilde5 = cache.tab

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
  u = @. uprev + γ*z₁
  b = dt.*f(tstep,u) .- z₁
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₁ = z₁ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. uprev + γ*z₁
    b = dt.*f(tstep,u) .- z₁
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
    z₁ = z₁ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ##### Step 2

  z₂ = α21.*z₁

  # initial step of Newton iteration
  iter = 1
  tstep = t + c2*dt
  tmp = @. uprev + a21*z₁
  u = @. tmp + γ*z₂
  b = dt.*f(tstep,u) .- z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₂
    b = dt.*f(tstep,u) .- z₂
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
    z₂ = z₂ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  z₃ = @. α31*z₁ + α32*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  tmp = @. uprev + a31*z₁ + a32*z₂
  u = @. tmp + γ*z₃
  b = dt.*f(tstep,u) .- z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₃
    b = dt.*f(tstep,u) .- z₃
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
    z₃ = z₃ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  z₄ = @. α41*z₁ + α43*z₃

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  tmp = @. uprev + a41*z₁ + a42*z₂ + a43*z₃
  u = @. tmp + γ*z₄
  b = dt.*f(tstep,u) .- z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₄
    b = dt.*f(tstep,u) .- z₄
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
    z₄ = z₄ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat2 for prediction
  z₅ = @. bhat1*z₁ + bhat2*z₂ + bhat3*z₃ + bhat4*z₄

  # initial step of Newton iteration
  iter = 1
  tstep = t+dt
  tmp = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  u = @. tmp + γ*z₅
  b = dt.*f(tstep,u) .- z₅
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₅ = z₅ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₅
    b = dt.*f(tstep,u) .- z₅
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
    z₅ = z₅ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. tmp + γ*z₅

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
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
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::Hairer4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,z₅,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,c2,c3,c4 = cache.tab
  @unpack α21,α31,α32,α41,α43 = cache.tab
  @unpack bhat1,bhat2,bhat3,bhat4,btilde1,btilde2,btilde3,btilde4,btilde5 = cache.tab

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,γdt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      mass_matrix = integrator.sol.prob.mass_matrix
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
      end
    else
      new_W = false
    end
  end

  # initial guess
  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
    @. z₁ = u - uprev
  elseif integrator.alg.extrapolant == :linear
    @. z₁ = dt*integrator.fsalfirst
  else
    z₁ .= zero(0)
  end

  ##### Step 1

  # initial step of Newton iteration
  iter = 1
  tstep = t + γdt
  @. u = uprev + γ*z₁
  f(tstep,u,k)
  @. b = dt*k - z₁
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₁ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + γ*z₁
    f(tstep,u,k)
    @. b = dt*k - z₁
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₁ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(tstep,u,k)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₂
    f(tstep,u,k)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₂ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(tstep,u,k)
  @. b = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₃ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₃
    f(tstep,u,k)
    @. b = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₃ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(tstep,u,k)
  @. b = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₄ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₄
    f(tstep,u,k)
    @. b = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₄ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(tstep,u,k)
  @. b = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₅ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₅
    f(tstep,u,k)
    @. b = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₅ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        integrator.alg.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z₅/dt
end

function initialize!(integrator, cache::Kvaerno4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end


@muladd function perform_step!(integrator, cache::Kvaerno4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,c3,c4 = cache.tab
  @unpack α21,α31,α32,α41,α42 = cache.tab
  @unpack btilde1,btilde2,btilde3,btilde4,btilde5 = cache.tab

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

  z₁ = dt.*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation choice
  z₂ = zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt
  tmp = @. uprev + γ*z₁
  u = @. tmp + γ*z₂
  b = dt.*f(tstep,u) .- z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₂
    b = dt.*f(tstep,u) .- z₂
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
    z₂ = z₂ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  z₃ = @. α31*z₁ + α32*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  tmp = @. uprev + a31*z₁ + a32*z₂
  u = @. tmp + γ*z₃
  b = dt.*f(tstep,u) .- z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₃
    b = dt.*f(tstep,u) .- z₃
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
    z₃ = z₃ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  z₄ = @. α41*z₁ + α42*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  tmp = @. uprev + a41*z₁ + a42*z₂ + a43*z₃
  u = @. tmp + γ*z₄
  b = dt.*f(tstep,u) .- z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₄
    b = dt.*f(tstep,u) .- z₄
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
    z₄ = z₄ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat2 for prediction
  z₅ = @. a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  u = @. tmp + γ*z₅
  b = dt.*f(tstep,u) .- z₅
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₅ = z₅ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₅
    b = dt.*f(tstep,u) .- z₅
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
    z₅ = z₅ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. tmp + γ*z₅

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = z₅./dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::Kvaerno4Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::Kvaerno4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,z₅,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,c3,c4 = cache.tab
  @unpack α21,α31,α32,α41,α42 = cache.tab
  @unpack btilde1,btilde2,btilde3,btilde4,btilde5 = cache.tab

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,γdt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      mass_matrix = integrator.sol.prob.mass_matrix
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##### Step 1

  @. z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Allow other choices here
  z₂ .= zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt
  @. tmp = uprev + γ*z₁
  @. u = tmp + γ*z₂
  f(tstep,u,k)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₂
    f(tstep,u,k)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₂ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(tstep,u,k)
  @. b = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₃ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₃
    f(tstep,u,k)
    @. b = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₃ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = α41*z₁ + α42*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  @. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃
  @. u = tmp + γ*z₄
  f(tstep,u,k)
  @. b = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₄ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₄
    f(tstep,u,k)
    @. b = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₄ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  # Use yhat prediction
  @. z₅ = a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. tmp = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  @. u = tmp + γ*z₅
  f(tstep,u,k)
  @. b = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₅ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₅
    f(tstep,u,k)
    @. b = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₅ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + γ*z₅

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    @. dz = btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        integrator.alg.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z₅/dt
end

function initialize!(integrator, cache::KenCarp4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::KenCarp4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c3,c4,c5 = cache.tab
  @unpack α31,α32,α41,α42,α51,α52,α53,α54,α61,α62,α63,α64,α65 = cache.tab
  @unpack btilde1,btilde3,btilde4,btilde5,btilde6 = cache.tab

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

  z₁ = dt.*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation choice
  z₂ = zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt
  tmp = @. uprev + γ*z₁
  u = @. tmp + γ*z₂
  b = dt.*f(tstep,u) .- z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₂
    b = dt.*f(tstep,u) .- z₂
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
    z₂ = z₂ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  z₃ = @. α31*z₁ + α32*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  tmp = @. uprev + a31*z₁ + a32*z₂
  u = @. tmp + γ*z₃
  b = dt.*f(tstep,u) .- z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₃
    b = dt.*f(tstep,u) .- z₃
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
    z₃ = z₃ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  z₄ = @. α41*z₁ + α42*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  tmp = @. uprev + a41*z₁ + a42*z₂ + a43*z₃
  u = @. tmp + γ*z₄
  b = dt.*f(tstep,u) .- z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₄
    b = dt.*f(tstep,u) .- z₄
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
    z₄ = z₄ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  z₅ = @. α51*z₁ + α52*z₂ + α53*z₃ + α54*z₄

  # initial step of Newton iteration
  iter = 1
  tstep = t + c5*dt
  tmp = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  u = @. tmp + γ*z₅
  b = dt.*f(tstep,u) .- z₅
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₅ = z₅ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₅
    b = dt.*f(tstep,u) .- z₅
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
    z₅ = z₅ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  z₆ = @. α61*z₁ + α62*z₂ + α63*z₃ + α64*z₄ + α65*z₅

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = @. uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅
  u = @. tmp + γ*z₆
  b = dt.*f(tstep,u) .- z₆
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₆ = z₆ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₆
    b = dt.*f(tstep,u) .- z₆
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
    z₆ = z₆ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. tmp + γ*z₆

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. btilde1*z₁ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅ + btilde6*z₆
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = z₆./dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::KenCarp4Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::KenCarp4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,z₅,z₆,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c3,c4,c5 = cache.tab
  @unpack α31,α32,α41,α42,α51,α52,α53,α54,α61,α62,α63,α64,α65 = cache.tab
  @unpack btilde1,btilde3,btilde4,btilde5,btilde6 = cache.tab

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,γdt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      mass_matrix = integrator.sol.prob.mass_matrix
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##### Step 1

  @. z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Allow other choices here
  z₂ .= zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt
  @. tmp = uprev + γ*z₁
  @. u = tmp + γ*z₂
  f(tstep,u,k)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₂
    f(tstep,u,k)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₂ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(tstep,u,k)
  @. b = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₃ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₃
    f(tstep,u,k)
    @. b = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₃ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = α41*z₁ + α42*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  @. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃
  @. u = tmp + γ*z₄
  f(tstep,u,k)
  @. b = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₄ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₄
    f(tstep,u,k)
    @. b = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₄ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  @. z₅ = α51*z₁ + α52*z₂ + α53*z₃ + α54*z₄

  # initial step of Newton iteration
  iter = 1
  tstep = t + c5*dt
  @. tmp = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  @. u = tmp + γ*z₅
  f(tstep,u,k)
  @. b = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₅ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₅
    f(tstep,u,k)
    @. b = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₅ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 6

  #@. z₆ = α61*z₁ + α62*z₂ + α63*z₃ + α64*z₄ + α65*z₅
  @tight_loop_macros for i in eachindex(u)
    @inbounds z₆[i] = α61*z₁[i] + α62*z₂[i] + α63*z₃[i] + α64*z₄[i] + α65*z₅[i]
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  @. tmp = uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅
  @. u = tmp + γ*z₆
  f(tstep,u,k)
  @. b = dt*k - z₆
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₆ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₆
    f(tstep,u,k)
    @. b = dt*k - z₆
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₆ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + γ*z₆

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    # @. dz = btilde1*z₁ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅ + btilde6*z₆
    @tight_loop_macros for i in eachindex(u)
      @inbounds dz[i] = btilde1*z₁[i] + btilde3*z₃[i] + btilde4*z₄[i] + btilde5*z₅[i] + btilde6*z₆[i]
    end
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        integrator.alg.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z₆/dt
end

function initialize!(integrator, cache::Kvaerno5ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Kvaerno5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,c3,c4,c5,c6 = cache.tab
  @unpack btilde1,btilde3,btilde4,btilde5,btilde6,btilde7 = cache.tab
  @unpack α31,α32,α41,α42,α43,α51,α52,α53,α61,α62,α63 = cache.tab

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

  z₁ = dt.*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation choice
  z₂ = zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt
  tmp = @. uprev + γ*z₁
  u = @. tmp + γ*z₂
  b = dt.*f(tstep,u) .- z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₂
    b = dt.*f(tstep,u) .- z₂
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
    z₂ = z₂ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  z₃ = @. α31*z₁ + α32*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  tmp = @. uprev + a31*z₁ + a32*z₂
  u = @. tmp + γ*z₃
  b = dt.*f(tstep,u) .- z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₃
    b = dt.*f(tstep,u) .- z₃
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
    z₃ = z₃ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  z₄ = @. α41*z₁ + α42*z₂ + α43*z₃

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  tmp = @. uprev + a41*z₁ + a42*z₂ + a43*z₃
  u = @. tmp + γ*z₄
  b = dt.*f(tstep,u) .- z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₄
    b = dt.*f(tstep,u) .- z₄
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
    z₄ = z₄ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  z₅ = @. α51*z₁ + α52*z₂ + α53*z₃

  # initial step of Newton iteration
  iter = 1
  tstep = t + c5*dt
  tmp = @. uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  u = @. tmp + γ*z₅
  b = dt.*f(tstep,u) .- z₅
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₅ = z₅ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₅
    b = dt.*f(tstep,u) .- z₅
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
    z₅ = z₅ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 6

  z₆ = @. α61*z₁ + α62*z₂ + α63*z₃

  # initial step of Newton iteration
  iter = 1
  tstep = t + c6*dt
  tmp = @. uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅
  u = @. tmp + γ*z₆
  b = dt.*f(tstep,u) .- z₆
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₆ = z₆ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₆
    b = dt.*f(tstep,u) .- z₆
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
    z₆ = z₆ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 7

  # Prediction from embedding
  z₇ = @. a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = @. uprev + a71*z₁ + a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆
  u = @. tmp + γ*z₇
  b = dt.*f(tstep,u) .- z₇
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₇ = z₇ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₇
    b = dt.*f(tstep,u) .- z₇
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
    z₇ = z₇ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. tmp + γ*z₇

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. btilde1*z₁ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅ + btilde6*z₆ + btilde7*z₇
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = z₇./dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::Kvaerno5Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::Kvaerno5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,z₅,z₆,z₇,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,c3,c4,c5,c6 = cache.tab
  @unpack btilde1,btilde3,btilde4,btilde5,btilde6,btilde7 = cache.tab
  @unpack α31,α32,α41,α42,α43,α51,α52,α53,α61,α62,α63 = cache.tab

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,γdt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      mass_matrix = integrator.sol.prob.mass_matrix
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##### Step 1

  @. z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Allow other choices here
  z₂ .= zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt
  @. tmp = uprev + γ*z₁
  @. u = tmp + γ*z₂
  f(tstep,u,k)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₂
    f(tstep,u,k)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₂ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(tstep,u,k)
  @. b = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₃ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₃
    f(tstep,u,k)
    @. b = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₃ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = α41*z₁ + α42*z₂ + α43*z₃

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  @. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃
  @. u = tmp + γ*z₄
  f(tstep,u,k)
  @. b = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₄ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₄
    f(tstep,u,k)
    @. b = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₄ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  @. z₅ = α51*z₁ + α52*z₂ + α53*z₃

  # initial step of Newton iteration
  iter = 1
  tstep = t + c5*dt
  @. tmp = uprev + a51*z₁ + a52*z₂ + a53*z₃ + a54*z₄
  @. u = tmp + γ*z₅
  f(tstep,u,k)
  @. b = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₅ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₅
    f(tstep,u,k)
    @. b = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₅ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 6

  @. z₆ = α61*z₁ + α62*z₂ + α63*z₃

  # inital step of Newton iteration
  iter = 1
  tstep = t + c6*dt
  @. tmp = uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅
  @. u = tmp + γ*z₆
  f(tstep,u,k)
  @. b = dt*k - z₆
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₆ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₆
    f(tstep,u,k)
    @. b = dt*k - z₆
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₆ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 7

  # Prediction is embedded method
  # @. z₇ = a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅ + γ*z₆
  @tight_loop_macros for i in eachindex(u)
    @inbounds z₇[i] = a61*z₁[i] + a63*z₃[i] + a64*z₄[i] + a65*z₅[i] + γ*z₆[i]
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  # @. tmp = uprev + a71*z₁ + a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆
  @tight_loop_macros for i in eachindex(u)
    @inbounds tmp[i] = uprev[i] + a71*z₁[i] + a73*z₃[i] + a74*z₄[i] + a75*z₅[i] + a76*z₆[i]
  end
  @. u = tmp + γ*z₇
  f(tstep,u,k)
  @. b = dt*k - z₇
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₇ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₇
    f(tstep,u,k)
    @. b = dt*k - z₇
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₇ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + γ*z₇

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    # @. dz = btilde1*z₁ + btilde3*z₃ + btilde4*z₄ + btilde5*z₅ + btilde6*z₆ + btilde7*z₇
    @tight_loop_macros for i in eachindex(u)
      @inbounds dz[i] = btilde1*z₁[i] + btilde3*z₃[i] + btilde4*z₄[i] + btilde5*z₅[i] + btilde6*z₆[i] + btilde7*z₇[i]
    end
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        integrator.alg.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z₇/dt
end

function initialize!(integrator, cache::KenCarp5ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::KenCarp5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a84,a85,a86,a87,c3,c4,c5,c6,c7 = cache.tab
  @unpack α31,α32,α41,α42,α51,α52,α61,α62,α71,α72,α73,α74,α75,α81,α82,α83,α84,α85 = cache.tab
  @unpack btilde1,btilde4,btilde5,btilde6,btilde7,btilde8 = cache.tab

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

  z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation choice
  z₂ = zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt
  tmp = @. uprev + γ*z₁
  u = @. tmp + γ*z₂
  b = dt.*f(tstep,u) .- z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ .+ dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₂
    b = dt.*f(tstep,u) .- z₂
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
    z₂ = z₂ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  z₃ = @. α31*z₁ + α32*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  tmp = @. uprev + a31*z₁ + a32*z₂
  u = @. tmp + γ*z₃
  b = dt.*f(tstep,u) .- z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₃
    b = dt.*f(tstep,u) .- z₃
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
    z₃ = z₃ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  z₄ = @. α41*z₁ + α42*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  tmp = @. uprev + a41*z₁ + a43*z₃
  u = @. tmp + γ*z₄
  b = dt.*f(tstep,u) .- z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₄
    b = dt.*f(tstep,u) .- z₄
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
    z₄ = z₄ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  z₅ = @. α51*z₁ + α52*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c5*dt
  tmp = @. uprev + a51*z₁ + a53*z₃ + a54*z₄
  u = @. tmp + γ*z₅
  b = dt.*f(tstep,u) .- z₅
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₅ = z₅ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₅
    b = dt.*f(tstep,u) .- z₅
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
    z₅ = z₅ + dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 6

  z₆ = @. α61*z₁ + α62*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c6*dt
  tmp = @. uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅
  u = @. tmp + γ*z₆
  b = dt.*f(tstep,u) .- z₆
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₆ = z₆ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₆
    b = dt.*f(tstep,u) .- z₆
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
    z₆ = z₆ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 7

  z₇ = @. α71*z₁ + α72*z₂ + α73*z₃ + α74*z₄ + α75*z₅

  # initial step of Newton iteration
  iter = 1
  tstep = t + c7*dt
  tmp = @. uprev + a71*z₁ +  a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆
  u = @. tmp + γ*z₇
  b = dt.*f(tstep,u) .- z₇
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₇ = z₇ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₇
    b = dt.*f(tstep,u) .- z₇
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
    z₇ = z₇ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 8

  # Prediction from embedding
  z₈ = @. α81*z₁ + α82*z₂ + α83*z₃ + α84*z₄ + α85*z₅

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  tmp = @. uprev + a81*z₁ + a84*z₄ + a85*z₅ + a86*z₆ + a87*z₇
  u = @. tmp + γ*z₈
  b = dt.*f(tstep,u) .- z₈
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₈ = z₈ .+ dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₈
    b = dt.*f(tstep,u) .- z₈
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
    z₈ = z₈ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. tmp + γ*z₈

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    # tmp = @. btilde1*z₁ + btilde4*z₄ + btilde5*z₅ + btilde6*z₆ + btilde7*z₇ + btilde8*z₈
    tmp = @. btilde1*z₁ + btilde4*z₄ + btilde5*z₅ + btilde6*z₆ + btilde7*z₇ + btilde8*z₈
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = z₈./dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::KenCarp5Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::KenCarp5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,z₅,z₆,z₇,z₈,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a84,a85,a86,a87,c3,c4,c5,c6,c7 = cache.tab
  @unpack α31,α32,α41,α42,α51,α52,α61,α62,α71,α72,α73,α74,α75,α81,α82,α83,α84,α85 = cache.tab
  @unpack btilde1,btilde4,btilde5,btilde6,btilde7,btilde8 = cache.tab

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  # calculate W
  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,γdt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, vec(uprev), vec(du1), integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(integrator.EEst)*oneunit(integrator.t))
      new_W = true
      mass_matrix = integrator.sol.prob.mass_matrix
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##### Step 1

  @. z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Allow other choices here
  z₂ .= zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt
  @. tmp = uprev + γ*z₁
  @. u = tmp + γ*z₂
  f(tstep,u,k)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(cache.ηold,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₂
    f(tstep,u,k)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₂ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
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
  f(tstep,u,k)
  @. b = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₃ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₃
    f(tstep,u,k)
    @. b = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₃ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  # Use constant z prediction
  @. z₄ = α41*z₁ + α42*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c4*dt
  @. tmp = uprev + a41*z₁ + a43*z₃
  @. u = tmp + γ*z₄
  f(tstep,u,k)
  @. b = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₄ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₄
    f(tstep,u,k)
    @. b = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₄ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 5

  @. z₅ = α51*z₁ + α52*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c5*dt
  @. tmp = uprev + a51*z₁ + a53*z₃ + a54*z₄
  @. u = tmp + γ*z₅
  f(tstep,u,k)
  @. b = dt*k - z₅
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₅ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₅
    f(tstep,u,k)
    @. b = dt*k - z₅
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₅ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 6

  @. z₆ = α61*z₁ + α62*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c6*dt
  @. tmp = uprev + a61*z₁ + a63*z₃ + a64*z₄ + a65*z₅
  @. u = tmp + γ*z₆
  f(tstep,u,k)
  @. b = dt*k - z₆
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₆ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₆
    f(tstep,u,k)
    @. b = dt*k - z₆
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₆ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 7

  #@. z₇ = α71*z₁ + α72*z₂ + α73*z₃ + α74*z₄ + α75*z₅
  @tight_loop_macros for i in eachindex(u)
    @inbounds z₇[i] = α71*z₁[i] + α72*z₂[i] + α73*z₃[i] + α74*z₄[i] + α75*z₅[i]
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + c7*dt
  #@. tmp = uprev + a71*z₁ + a73*z₃ + a74*z₄ + a75*z₅ + a76*z₆
  @tight_loop_macros for i in eachindex(u)
    @inbounds tmp[i] = uprev[i] + a71*z₁[i] + a73*z₃[i] + a74*z₄[i] + a75*z₅[i] + a76*z₆[i]
  end
  @. u = tmp + γ*z₇
  f(tstep,u,k)
  @. b = dt*k - z₇
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₇ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₇
    f(tstep,u,k)
    @. b = dt*k - z₇
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₇ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 8

  #@. z₈ = α81*z₁ + α82*z₂ + α83*z₃ + α84*z₄ + α85*z₅
  @tight_loop_macros for i in eachindex(u)
    @inbounds z₈[i] = α81*z₁[i] + α82*z₂[i] + α83*z₃[i] + α84*z₄[i] + α85*z₅[i]
  end

  # initial step of Newton iteration
  iter = 1
  tstep = t + dt
  #@. u = uprev + a81*z₁ + a84*z₄ + a85*z₅ + a86*z₆ + a87*z₇
  @tight_loop_macros for i in eachindex(u)
    @inbounds tmp[i] = uprev[i] + a81*z₁[i] + a84*z₄[i] + a85*z₅[i] + a86*z₆[i] + a87*z₇[i]
  end
  @. u = tmp + γ*z₈
  f(tstep,u,k)
  @. b = dt*k - z₈
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₈ .+= dz

  η = max(η,eps(first(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₈
    f(tstep,u,k)
    @. b = dt*k - z₈
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(b),false)
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
    z₈ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + γ*z₈

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    # @. dz = btilde1*z₁ + btilde4*z₄ + btilde5*z₅ + btilde6*z₆ + btilde7*z₇ + btilde8*z₈
    @tight_loop_macros for i in eachindex(u)
      @inbounds dz[i] = btilde1*z₁[i] + btilde4*z₄[i] + btilde5*z₅[i] + btilde6*z₆[i] + btilde7*z₇[i] + btilde8*z₈[i]
    end
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        integrator.alg.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z₈/dt
end
