function initialize!(integrator, cache::ABDF2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ABDF2ConstantCache, repeat_step=false)
  @unpack t,f,p = integrator
  @unpack uf,κ,tol,dtₙ₋₁ = cache
  alg = unwrap_alg(integrator, true)
  dtₙ, uₙ, uₙ₋₁, uₙ₋₂ = integrator.dt, integrator.u, integrator.uprev, integrator.uprev2

  if integrator.iter == 1 && !integrator.u_modified
    cache.dtₙ₋₁ = dtₙ
    perform_step!(integrator, cache.eulercache, repeat_step)
    cache.fsalfirstprev = integrator.fsalfirst
    return
  end

  # precalculations
  fₙ₋₁ = integrator.fsalfirst
  κtol = κ*tol
  ρ = dtₙ/dtₙ₋₁
  d = 2/3
  ddt = d*dtₙ
  dtmp = ρ^2/3
  d1 = 1+dtmp
  d2 = -dtmp
  d3 = -(ρ-1)/3

  # calculate W
  W = calc_W!(integrator, cache, ddt, repeat_step)

  zₙ₋₁ = dtₙ*fₙ₋₁
  # initial guess
  if alg.extrapolant == :linear
    z = dtₙ*fₙ₋₁
  else # :constant
    z = zero(uₙ)
  end

  tmp = d1*uₙ₋₁ + d2*uₙ₋₂ + d3*zₙ₋₁
  z,η,iter,fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, d, 1, Val{:newton})
  fail_convergence && return

  uₙ = tmp + d*z
  integrator.fsallast = f(uₙ,p,t+dtₙ)

  if integrator.opts.adaptive
    tmp = integrator.fsallast - (1+dtₙ/dtₙ₋₁)*integrator.fsalfirst + (dtₙ/dtₙ₋₁)*cache.fsalfirstprev
    est = (dtₙ₋₁+dtₙ)/6 * tmp
    atmp = calculate_residuals(est, uₙ₋₁, uₙ, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  ################################### Finalize

  cache.dtₙ₋₁ = dtₙ
  cache.ηold = η
  cache.newton_iters = iter
  if integrator.EEst < one(integrator.EEst)
    cache.fsalfirstprev = integrator.fsalfirst
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = uₙ
  return
end

function initialize!(integrator, cache::ABDF2Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::ABDF2Cache, repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack dz,z,k,b,J,W,tmp,atmp,κ,tol,dtₙ₋₁,zₙ₋₁ = cache
  alg = unwrap_alg(integrator, true)
  uₙ,uₙ₋₁,uₙ₋₂,dtₙ = integrator.u,integrator.uprev,integrator.uprev2,integrator.dt

  if integrator.iter == 1 && !integrator.u_modified
    cache.dtₙ₋₁ = dtₙ
    perform_step!(integrator, cache.eulercache, repeat_step)
    cache.fsalfirstprev .= integrator.fsalfirst
    return
  end

  # precalculations
  fₙ₋₁ = integrator.fsalfirst
  κtol = κ*tol
  ρ = dtₙ/dtₙ₋₁
  d = 2/3
  ddt = d*dtₙ
  dtmp = ρ^2/3
  d1 = 1+dtmp
  d2 = -dtmp
  d3 = -(ρ-1)/3

  new_W = calc_W!(integrator, cache, ddt, repeat_step)

  # initial guess
  @. zₙ₋₁ = dtₙ*fₙ₋₁
  if alg.extrapolant == :linear
    @. z = dtₙ*fₙ₋₁
  else # :constant
    fill!(z, zero(eltype(z)))
  end

  @. tmp = d1*uₙ₋₁ + d2*uₙ₋₂ + d3*zₙ₋₁
  z,η,iter,fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, d, 1, Val{:newton}, new_W)
  fail_convergence && return

  @. uₙ = tmp + d*z

  f(integrator.fsallast, uₙ, p, t+dtₙ)
  if integrator.opts.adaptive
    btilde0 = (dtₙ₋₁+dtₙ)/6
    btilde1 = 1+dtₙ/dtₙ₋₁
    btilde2 = dtₙ/dtₙ₋₁
    @. tmp = btilde0*(integrator.fsallast - btilde1*integrator.fsalfirst + btilde2*cache.fsalfirstprev)
    calculate_residuals!(atmp, tmp, uₙ₋₁, uₙ, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter
  cache.dtₙ₋₁ = dtₙ
  if integrator.EEst < one(integrator.EEst)
    @. cache.fsalfirstprev = integrator.fsalfirst
  end
  return
end

# QNDF1

function initialize!(integrator, cache::QNDF1ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::QNDF1ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uprev2,D,D2,R,U,dtₙ₋₁ = cache
  cnt = integrator.iter
  k = 1
  if cnt == 1
    κ = zero(integrator.alg.kappa)
  else
    κ = integrator.alg.kappa
    ρ = dt/dtₙ₋₁
    D[1] = uprev - uprev2   # backward diff
    if ρ != 1
      R!(k,ρ,cache)
      D[1] = D[1] * (R[1] * U[1])
    end
  end

  # precalculations
  γ₁ = 1//1
  γ = inv((1-κ)*γ₁)

  u₀ = uprev + D[1]
  ϕ = γ * (γ₁*D[1])
  tmp = u₀ - ϕ

  γdt = γ*dt
  W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  z = dt*integrator.fsalfirst

  z, η, iter, fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, γ, 1, Val{:newton})
  fail_convergence && return
  u = tmp + γ*z

  if integrator.opts.adaptive && integrator.success_iter > 0
    D2[1] = u - uprev
    D2[2] = D2[1] - D[1]
    utilde = (κ*γ₁ + inv(k+1)) * D2[2]
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
    if integrator.EEst > one(integrator.EEst)
      return
    end
  else
    integrator.EEst = one(integrator.EEst)
  end
  cache.dtₙ₋₁ = dt
  cache.uprev2 = uprev
  cache.ηold = η
  cache.newton_iters = iter
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::QNDF1Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

function perform_step!(integrator,cache::QNDF1Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uprev2,D,D2,R,U,dtₙ₋₁,tmp,z,W,utilde,atmp = cache
  cnt = integrator.iter
  k = 1
  if cnt == 1
    κ = zero(integrator.alg.kappa)
  else
    κ = integrator.alg.kappa
    ρ = dt/dtₙ₋₁
    @. D[1] = uprev - uprev2 # backward diff
    if ρ != 1
      R!(k,ρ,cache)
      @. D[1] = D[1] * (R[1] * U[1])
    end
  end

  # precalculations
  γ₁ = 1//1
  γ = inv((1-κ)*γ₁)
  @. tmp = uprev + D[1] - γ * (γ₁*D[1])

  γdt = γ*dt
  new_W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  @. z = dt*integrator.fsalfirst

  z, η, iter, fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, γ, 1, Val{:newton}, new_W)
  fail_convergence && return
  @. u = tmp + γ*z

  if integrator.opts.adaptive && integrator.success_iter > 0
    @. D2[1] = u - uprev
    @. D2[2] = D2[1] - D[1]
    @. utilde = (κ*γ₁ + inv(k+1)) * D2[2]
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
    if integrator.EEst > one(integrator.EEst)
      return
    end
  else
    integrator.EEst = one(integrator.EEst)
  end
  cache.dtₙ₋₁ = dt
  cache.uprev2 .= uprev
  cache.ηold = η
  cache.newton_iters = iter
  f(integrator.fsallast, u, p, t+dt)
end

function initialize!(integrator, cache::QNDF2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::QNDF2ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uprev2,uprev3,dtₙ₋₁,dtₙ₋₂,D,D2,R,U = cache
  cnt = integrator.iter
  k = 2
  if cnt == 1 || cnt == 2
    κ = zero(integrator.alg.kappa)
    γ₁ = 1//1
    γ₂ = 1//1
  elseif dtₙ₋₁ != dtₙ₋₂
    κ = integrator.alg.kappa
    γ₁ = 1//1
    γ₂ = 1//1 + 1//2
    ρ₁ = dt/dtₙ₋₁
    ρ₂ = dt/dtₙ₋₂
    D[1] = uprev - uprev2
    D[1] = D[1] * ρ₁
    D[2] = D[1] - ((uprev2 - uprev3) * ρ₂)
  else
    κ = integrator.alg.kappa
    γ₁ = 1//1
    γ₂ = 1//1 + 1//2
    ρ = dt/dtₙ₋₁
    # backward diff
    D[1] = uprev - uprev2
    D[2] = D[1] - (uprev2 - uprev3)
    if ρ != 1
      R!(k,ρ,cache)
      R .= R * U
      D[1] = D[1] * R[1,1] + D[2] * R[2,1]
      D[2] = D[1] * R[1,2] + D[2] * R[2,2]
    end
  end

  # precalculations
  γ = inv((1-κ)*γ₂)
  u₀ = uprev + D[1] + D[2]
  ϕ = γ * (γ₁*D[1] + γ₂*D[2])
  tmp = u₀ - ϕ

  γdt = γ*dt
  W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  z = dt*integrator.fsalfirst

  z, η, iter, fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, γ, 1, Val{:newton})
  fail_convergence && return
  u = tmp + γ*z

  if integrator.opts.adaptive
    if integrator.success_iter == 0
      integrator.EEst = one(integrator.EEst)
    elseif integrator.success_iter == 1
      utilde = (u - uprev) - ((uprev - uprev2) * dt/dtₙ₋₁)
      atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp)
    else
      D2[1] = u - uprev
      D2[2] = D2[1] - D[1]
      D2[3] = D2[2] - D[2]
      utilde = (κ*γ₂ + inv(k+1)) * D2[3]
      atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp)
    end
  end
  if integrator.EEst > one(integrator.EEst)
    return
  end

  cache.uprev3 = uprev2
  cache.uprev2 = uprev
  cache.dtₙ₋₂ = dtₙ₋₁
  cache.dtₙ₋₁ = dt
  cache.ηold = η
  cache.newton_iters = iter
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::QNDF2Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

function perform_step!(integrator,cache::QNDF2Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uprev2,uprev3,dtₙ₋₁,dtₙ₋₂,D,D2,R,U,tmp,utilde,atmp,W,z = cache
  cnt = integrator.iter
  k = 2
  if cnt == 1 || cnt == 2
    κ = zero(integrator.alg.kappa)
    γ₁ = 1//1
    γ₂ = 1//1
  elseif dtₙ₋₁ != dtₙ₋₂
    κ = integrator.alg.kappa
    γ₁ = 1//1
    γ₂ = 1//1 + 1//2
    ρ₁ = dt/dtₙ₋₁
    ρ₂ = dt/dtₙ₋₂
    @. D[1] = uprev - uprev2
    @. D[1] = D[1] * ρ₁
    @. D[2] = D[1] - ((uprev2 - uprev3) * ρ₂)
  else
    κ = integrator.alg.kappa
    γ₁ = 1//1
    γ₂ = 1//1 + 1//2
    ρ = dt/dtₙ₋₁
    # backward diff
    @. D[1] = uprev - uprev2
    @. D[2] = D[1] - (uprev2 - uprev3)
    if ρ != 1
      R!(k,ρ,cache)
      R .= R * U
      @. D[1] = D[1] * R[1,1] + D[2] * R[2,1]
      @. D[2] = D[1] * R[1,2] + D[2] * R[2,2]
    end
  end

  # precalculations
  γ = inv((1-κ)*γ₂)
  @. tmp = uprev + D[1] + D[2] - γ * (γ₁*D[1] + γ₂*D[2])

  γdt = γ*dt
  new_W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  @. z = dt*integrator.fsalfirst

  z, η, iter, fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, γ, 1, Val{:newton}, new_W)
  fail_convergence && return
  @. u = tmp + γ*z

  if integrator.opts.adaptive
    if integrator.success_iter == 0
      integrator.EEst = one(integrator.EEst)
    elseif integrator.success_iter == 1
      @. utilde = (u - uprev) - ((uprev - uprev2) * dt/dtₙ₋₁)
      calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp)
    else
      @. D2[1] = u - uprev
      @. D2[2] = D2[1] - D[1]
      @. D2[3] = D2[2] - D[2]
      @. utilde = (κ*γ₂ + inv(k+1)) * D2[3]
      calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp)
    end
  end
  if integrator.EEst > one(integrator.EEst)
    return
  end

  cache.uprev3 .= uprev2
  cache.uprev2 .= uprev
  cache.dtₙ₋₂ = dtₙ₋₁
  cache.dtₙ₋₁ = dt
  cache.ηold = η
  cache.newton_iters = iter
  f(integrator.fsallast, u, p, t+dt)
end

function initialize!(integrator, cache::QNDFConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::QNDFConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack udiff,dts,k,max_order,D,D2,R,U = cache
  cnt = integrator.iter
  κ = integrator.alg.kappa[k]
  γ = inv((1-κ)*γₖ[k])
  flag = true
  for i in 2:k
    if dts[i] != dts[1]
      flag = false
      break
    end
  end

  if cnt > 2
    if flag
      ρ = dt/dts[1]
      # backward diff
      backward_diff(udiff,D,D2,k)
      if ρ != 1
        U!(k,U)
        R!(k,ρ,cache)
        R .= R * U
        D .= D * R
      end
    else
      for i = 1:k
        udiff[i] *= dt/dts[i]
      end
      backward_diff(udiff,D,D2,k)
    end
  else
    γ = 1//1
  end
  # precalculations
  u₀ = uprev + sum(D)
  ϕ = zero(γ)
  for i in 1:k
    ϕ += γₖ[k]*D[i]
  end
  ϕ *= γ
  tmp = u₀ - ϕ

  γdt = γ*dt
  W = calc_W!(integrator, cache, γdt, repeat_step)

  # initial guess
  z = dt*integrator.fsalfirst

  z, η, iter, fail_convergence = diffeq_nlsolve!(integrator, cache, W, z, tmp, γ, 1, Val{:newton})
  fail_convergence && return
  u = tmp + γ*z

  if integrator.opts.adaptive
    if integrator.success_iter == 0
      integrator.EEst = one(integrator.EEst)
    elseif integrator.success_iter == 1
      utilde = (u - uprev) - (udiff[1] * dt/dtₙ₋₁)
      atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp)
    else
      δ = u - uprev
      for i = 1:k
        δ -= δ - D[i]
      end
      utilde = (κ*γₖ[k] + inv(k+1)) * δ
      atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp)
    end
    if integrator.EEst > one(integrator.EEst)
      return
    end

    if cnt <=  4 || k < 3
      cache.k = min(k+1,3)
      if cnt == 1
        cache.k = 1
      end
    else
      utildem1 = (κ*γₖ[k-1] + inv(k)) * D[k]
      utildem2 = (κ*γₖ[k-2] + inv(k-1)) * D[k-1]
      backward_diff(udiff,D,D2,k+1)
      δ = u - uprev
      for i = 1:(k+1)
        δ -= δ - D[i]
      end
      utildep1 = (κ*γₖ[k+1] + inv(k+2)) * δ
      atmpm2 = calculate_residuals(utildem2, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      atmpm1 = calculate_residuals(utildem1, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      atmpp1 = calculate_residuals(utildep1, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      errm2 = integrator.opts.internalnorm(atmpm2)
      errm1 = integrator.opts.internalnorm(atmpm1)
      errp1 = integrator.opts.internalnorm(atmpp1)
      if max(errm2,errm1) <= integrator.EEst
        cache.k = k - 1
      elseif errp1 < integrator.EEst
        cache.k = min(k+1,max_order)
        integrator.EEst = one(integrator.EEst)   # for keeping the stepsize constant in the next step
      end # if
    end # step <= 4
  end # integrator.opts.adaptive

  for i = 2:5
    dts[i] = dts[i-1]
    udiff[i] = udiff[i-1]
  end
  dts[1] = dt
  udiff[1] = u - uprev

  cache.ηold = η
  cache.newton_iters = iter
  integrator.fsallast = f(u, p, t+dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end
