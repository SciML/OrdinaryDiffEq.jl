function initialize!(integrator, cache::ImplicitEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ImplicitEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache

  κtol = κ*tol

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    u = current_extrapolant(t+dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    u = @. uprev + integrator.fsalfirst*dt
  else # :constant
    u = uprev
  end

  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*J
  end

  z = u .- uprev
  iter = 1
  u = uprev .+ z
  b = dt.*f(t+dt,u) .- z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z .+ dz

  η = max(cache.ηold,eps(first(u)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = uprev .+ z
    b = dt.*f(t+dt,u) .- z
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

  cache.ηold = η
  cache.newton_iters = iter
  u = uprev .+ z

  if integrator.opts.adaptive && integrator.success_iter > 0
    # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
    # use 2nd divided difference (DD) a la SPICE and Shampine
    uprev2 = integrator.uprev2
    tprev = integrator.tprev

    # TODO: change recursive calculation to better algorithm
    diff10 = t-tprev
    diff11 = dt
    diff20 = t+dt-tprev
    c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
    r = c*dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

    dEst = @. r*abs(((u-uprev)/diff11-(uprev-uprev2)/diff10)/diff20)

    atmp = calculate_residuals(dEst, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol)
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
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::ImplicitEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z,k,J,W,jac_config,κ,tol = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  κtol = κ*tol

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    @. u = uprev + integrator.fsalfirst*dt
  else # :constant
    copy!(u,uprev)
  end

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
        if alg_autodiff(integrator.alg)
          ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
        else
          Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
        end
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps())
      new_W = true
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-dt*J[i,j]
      end
    else
      new_W = false
    end
  end

  @. z = u - uprev
  iter = 1
  @. u = uprev + z
  f(t+dt,u,k)
  scale!(k,dt)
  k .-= z
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz

  η = max(cache.ηold,eps(first(u)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + z
    f(t+dt,u,k)
    scale!(k,dt)
    k .-= z
    if has_invW(f)
      A_mul_B!(dz,W,k) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(k),false)
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
  @. u = uprev + z

  if integrator.opts.adaptive && integrator.success_iter > 0
    # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
    # use 2nd divided difference (DD) a la SPICE and Shampine
    uprev2 = integrator.uprev2
    tprev = integrator.tprev

    # TODO: change recursive calculation to better algorithm
    diff10 = t-tprev
    diff11 = dt
    diff20 = t+dt-tprev
    c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
    r = c*dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

    @. k = r*abs(((u-uprev)/diff11-(uprev-uprev2)/diff10)/diff20)

    # does not work with units - additional unitless array required!
    calculate_residuals!(k, k, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(k)
  else
    integrator.EEst = 1
  end

  f(t+dt,u,integrator.fsallast)
end

function initialize!(integrator, cache::TrapezoidConstantCache)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::TrapezoidConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache
  dto2 = dt/2

  κtol = κ*tol

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    u = current_extrapolant(t+dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    u = @. uprev + integrator.fsalfirst*dt
  else # :constant
    u = uprev
  end

  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dto2*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dto2*J
  end

  z = u .- uprev
  iter = 1
  u = uprev .+ z
  b = dto2.*f(t+dto2,u) .- z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z .+ dz

  η = max(cache.ηold,eps(first(u)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = uprev .+ z
    b = dto2.*f(t+dto2,u) .- z
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

  cache.ηold = η
  cache.newton_iters = iter
  u = @. uprev + 2*z

  if integrator.opts.adaptive
    if integrator.iter > 2
      # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
      # use 3rd divided differences (DD) a la SPICE and Shampine
      uprev2 = integrator.uprev2
      tprev = integrator.tprev
      uprev3 = cache.uprev3
      tprev2 = cache.tprev2

      # TODO: change recursive calculation to better algorithm
      diff10 = tprev-tprev2
      diff11 = t-tprev
      diff12 = dt
      diff20 = t-tprev2
      diff21 = t+dt-tprev
      diff30 = t+dt-tprev2
      c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
      r = c*dt^3/2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

      DD31 = @. ((u-uprev)/diff12 - (uprev-uprev2)/diff11)/diff21
      DD30 = @. ((uprev-uprev2)/diff11 - (uprev2-uprev3)/diff10)/diff20
      dEst = @. r*abs((DD31-DD30)/diff30) # scaled 3rd DD

      atmp = calculate_residuals(dEst, uprev, u, integrator.opts.abstol,
                                 integrator.opts.reltol)
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
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::TrapezoidCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z,k,J,W,jac_config,κ,tol = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  κtol = κ*tol

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    @. u = uprev + integrator.fsalfirst*dt
  else # :constant
    copy!(u,uprev)
  end

  dto2 = dt/2

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
        if alg_autodiff(integrator.alg)
          ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
        else
          Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
        end
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps())
      new_W = true
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-dto2*J[i,j]
      end
    else
      new_W = false
    end
  end

  @. z = u - uprev
  iter = 1
  @. u = uprev + z
  f(t+dto2,u,k)
  scale!(k,dto2)
  k .-= z
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz

  η = max(cache.ηold,eps(first(u)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol


  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = uprev + z
    f(t+dto2,u,k)
    scale!(k,dto2)
    k .-= z
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(k),false)
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
  @. u = uprev + 2*z

  if integrator.opts.adaptive
    if integrator.iter > 2
      # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
      # use 3rd divided difference (DD) a la SPICE and Shampine
      uprev2 = integrator.uprev2
      tprev = integrator.tprev
      uprev3 = cache.uprev3
      tprev2 = cache.tprev2

      # TODO: change recursive calculation to better algorithm
      diff10 = tprev-tprev2
      diff11 = t-tprev
      diff12 = dt
      diff20 = t-tprev2
      diff21 = t+dt-tprev
      diff30 = t+dt-tprev2
      c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
      r = c*dt^3/2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

      @. z = ((u-uprev)/diff12 - (uprev-uprev2)/diff11)/diff21 # DD31
      @. k = ((uprev-uprev2)/diff11 - (uprev2-uprev3)/diff10)/diff20 # DD30
      @. k = r*abs((z-k)/diff30) # scaled 3rd DD

      # does not work with units - additional unitless array required!
      calculate_residuals!(k, k, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
      integrator.EEst = integrator.opts.internalnorm(k)
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
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::TRBDF2ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache

  # Precalculations
  # -> move constants to cache?
  γ = 2 - sqrt(2)
  ω = sqrt(2)/4
  d = γ/2

  a1 = -sqrt(2)/2
  a2 = 1 + sqrt(2)/2

  btilde1 = (1-sqrt(2))/3
  btilde2 = 1/3
  btilde3 = (sqrt(2)-2)/3

  κtol = κ*tol

  γdt = γ*dt
  ddt = d*dt

  uf.t = t
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - ddt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - ddt*J
  end

  zprev = dt.*integrator.fsalfirst

  ##### Solve Trapezoid Step

  # TODO: Add extrapolant initial guess
  zᵧ = zprev # might change hence equality not used in calculations below

  iter = 1
  tmp = @. uprev + d*zprev # reduces calculations in Newton iteration
  uᵧ = @. tmp + d*zᵧ
  b = dt.*f(t+γdt,uᵧ) .- zᵧ
  Δzᵧ = W\b
  ndz = integrator.opts.internalnorm(Δzᵧ)
  zᵧ = zᵧ .+ Δzᵧ

  η = max(cache.ηold,eps(first(u)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    uᵧ = @. tmp + d*zᵧ
    b = dt.*f(t+γdt,uᵧ) .- zᵧ
    Δzᵧ = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(Δzᵧ)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    zᵧ = zᵧ .+ Δzᵧ
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve BDF2 Step

  ### Initial Guess From Shampine

  #=
  Original assignments
    uᵧ = @. (uprev + d*zprev) + d*zᵧ
    z = @. (1.5 + sqrt(2))*zprev + (2.5 + 2sqrt(2))*zᵧ - (6 + 4.5sqrt(2))*(uᵧ - uprev)
  with d = (2 - sqrt(2))/2 can be simplified to
    z = @. -sqrt(2)/2*zprev + (1 + sqrt(2)/2)*zᵧ
  =#
  z = @. a1*zprev + a2*zᵧ

  iter = 1
  tmp = @. uprev + ω*zprev + ω*zᵧ # reduces calculations in Newton iteration
  u = @. tmp + d*z
  b = dt.*f(t+dt,u) .- z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z .+ dz

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + d*z
    b = dt.*f(t+dt,u) .- z
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

  if integrator.opts.adaptive
    tmp = @. btilde1*zprev + btilde2*zᵧ + btilde3*z
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
                               integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = z./dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.ηold = η
  cache.newton_iters = iter
  integrator.u = u
end

function initialize!(integrator, cache::TRBDF2Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::TRBDF2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,uᵧ,Δzᵧ,Δz,zprev,zᵧ,z,k,J,W,jac_config,tmp,κ,tol = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  # Precalculations
  # -> move constants to cache?
  γ = 2 - sqrt(2)
  ω = sqrt(2)/4
  d = γ/2

  a1 = -sqrt(2)/2
  a2 = 1 + sqrt(2)/2

  btilde1 = (1-sqrt(2))/3
  btilde2 = 1/3
  btilde3 = (sqrt(2)-2)/3

  κtol = κ*tol

  γdt = γ*dt
  ddt = d*dt

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
        if alg_autodiff(integrator.alg)
          ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
        else
          Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
        end
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps())
      new_W = true
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-ddt*J[i,j]
      end
    else
      new_W = false
    end
  end

  @. zprev = dt*integrator.fsalfirst

  ##### Solve Trapezoid Step

  # TODO: Add extrapolant initial guess
  @. zᵧ = zprev # might change hence equality not used in calculations below

  iter = 1
  @. tmp = uprev + d*zprev # reduces calculations in Newton iteration
  @. uᵧ = tmp + d*zᵧ
  f(t+γdt,uᵧ,k)
  @. k = dt*k - zᵧ
  if has_invW(f)
    A_mul_B!(vec(Δzᵧ),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(Δzᵧ),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(Δzᵧ)
  zᵧ .+= Δzᵧ

  η = max(cache.ηold,eps(first(u)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. uᵧ = tmp + d*zᵧ
    f(t+γdt,uᵧ,k)
    @. k = dt*k - zᵧ
    if has_invW(f)
      A_mul_B!(vec(Δzᵧ),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(Δzᵧ),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(Δzᵧ)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    zᵧ .+= Δzᵧ
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve BDF2 Step

  ### Initial Guess From Shampine

  #=
  Original assignments
    @. uᵧ = (uprev + d*zprev) + d*zᵧ
    @. z = (1.5 + sqrt(2))*zprev + (2.5 + 2sqrt(2))*zᵧ - (6 + 4.5sqrt(2))*(uᵧ - uprev)
  with d = (2 - sqrt(2))/2 can be simplified to
    @. z = -sqrt(2)/2*zprev + (1 + sqrt(2)/2)*zᵧ
  =#
  @. z = a1*zprev + a2*zᵧ

  iter = 1
  @. tmp = uprev + ω*zprev + ω*zᵧ # reduces calculations in Newton iteration
  @. u = tmp + d*z
  f(t+dt,u,k)
  @. k = dt*k - z
  if has_invW(f)
    A_mul_B!(vec(Δz),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(Δz),W,vec(k),false)
  end
  ndz = integrator.opts.internalnorm(Δz)
  z .+= Δz

  η = max(η,eps(first(u)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + d*z
    f(t+dt,u,k)
    @. k = dt*k - z
    if has_invW(f)
      A_mul_B!(vec(Δz),W,vec(k)) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(Δz),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(Δz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z .+= Δz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  @. u = tmp + d*z

  ################################### Finalize

  if integrator.opts.adaptive
    @. tmp = btilde1*zprev + btilde2*zᵧ + btilde3*z
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(k),W,vec(tmp))
      else
        integrator.alg.linsolve(vec(k),W,vec(tmp),false)
      end
    else
      k .= tmp
    end
    # does not work with units - additional unitless array required!
    calculate_residuals!(tmp, k, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

  @. integrator.fsallast = z/dt
  cache.ηold = η
  cache.newton_iters = iter
end
