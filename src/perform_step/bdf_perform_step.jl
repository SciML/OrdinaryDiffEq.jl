function initialize!(integrator, cache::ABDF2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ABDF2ConstantCache, repeat_step=false)
  @unpack t,f,p = integrator
  @unpack dtₙ₋₁,nlsolver = cache
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
  ρ = dtₙ/dtₙ₋₁
  d = 2//3
  ddt = d*dtₙ
  dtmp = 1//3*ρ^2
  d1 = 1+dtmp
  d2 = -dtmp
  d3 = -(ρ-1)*1//3

  # calculate W
  markfirststage!(nlsolver)

  zₙ₋₁ = dtₙ*fₙ₋₁
  # initial guess
  if alg.extrapolant == :linear
    z = dtₙ*fₙ₋₁
  else # :constant
    z = zero(uₙ)
  end
  nlsolver.z = z

  nlsolver.tmp = d1*uₙ₋₁ + d2*uₙ₋₂ + d3*zₙ₋₁
  nlsolver.γ = d
  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return

  uₙ = nlsolver.tmp + d*z
  integrator.fsallast = f(uₙ,p,t+dtₙ)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    tmp = integrator.fsallast - (1+dtₙ/dtₙ₋₁)*integrator.fsalfirst + (dtₙ/dtₙ₋₁)*cache.fsalfirstprev
    est = (dtₙ₋₁+dtₙ)/6 * tmp
    atmp = calculate_residuals(est, uₙ₋₁, uₙ, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  ################################### Finalize

  if integrator.EEst < one(integrator.EEst)
    cache.fsalfirstprev = integrator.fsalfirst
    cache.dtₙ₋₁ = dtₙ
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = uₙ
  return
end

function initialize!(integrator, cache::ABDF2Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::ABDF2Cache, repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack atmp,dtₙ₋₁,zₙ₋₁,nlsolver = cache
  @unpack z,tmp = nlsolver
  alg = unwrap_alg(integrator, true)
  uₙ,uₙ₋₁,uₙ₋₂,dtₙ = integrator.u,integrator.uprev,integrator.uprev2,integrator.dt

  if integrator.iter == 1 && !integrator.u_modified
    cache.dtₙ₋₁ = dtₙ
    perform_step!(integrator, cache.eulercache, repeat_step)
    cache.fsalfirstprev .= integrator.fsalfirst
    nlsolver.tmp = tmp
    return
  end

  # precalculations
  fₙ₋₁ = integrator.fsalfirst
  ρ = dtₙ/dtₙ₋₁
  d = 2//3
  ddt = d*dtₙ
  dtmp = 1//3*ρ^2
  d1 = 1+dtmp
  d2 = -dtmp
  d3 = -(ρ-1)*1//3

  markfirststage!(nlsolver)

  @.. zₙ₋₁ = dtₙ*fₙ₋₁
  # initial guess
  if alg.extrapolant == :linear
    @.. z = dtₙ*fₙ₋₁
  else # :constant
    z .= zero(eltype(z))
  end

  @.. tmp = d1*uₙ₋₁ + d2*uₙ₋₂ + d3*zₙ₋₁
  nlsolver.γ = d
  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return

  @.. uₙ = tmp + d*z

  f(integrator.fsallast, uₙ, p, t+dtₙ)
  integrator.destats.nf += 1
  if integrator.opts.adaptive
    btilde0 = (dtₙ₋₁+dtₙ)*1//6
    btilde1 = 1+dtₙ/dtₙ₋₁
    btilde2 = dtₙ/dtₙ₋₁
    @.. tmp = btilde0*(integrator.fsallast - btilde1*integrator.fsalfirst + btilde2*cache.fsalfirstprev)
    calculate_residuals!(atmp, tmp, uₙ₋₁, uₙ, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  ################################### Finalize

  if integrator.EEst < one(integrator.EEst)
    @.. cache.fsalfirstprev = integrator.fsalfirst
    cache.dtₙ₋₁ = dtₙ
  end
  return
end

# SBDF

function initialize!(integrator, cache::SBDFConstantCache)
  @unpack uprev, p, t = integrator
  @unpack f1, f2 = integrator.f
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  cache.du₁ = f1(uprev,p,t)
  cache.du₂ = f2(uprev,p,t)
  integrator.destats.nf += 1
  integrator.destats.nf2 += 1
  integrator.fsalfirst = cache.du₁ + cache.du₂

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::SBDFConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p,alg = integrator
  @unpack uprev2,uprev3,uprev4,du₁,du₂,k₁,k₂,k₃,nlsolver = cache
  @unpack f1, f2 = integrator.f
  cnt = cache.cnt = min(alg.order, integrator.iter+1)
  integrator.iter == 1 && !integrator.u_modified && ( cnt = cache.cnt = 1 )
  nlsolver.γ = γ = inv(γₖ[cnt])
  if cnt == 1
    tmp = uprev + dt*du₂
  elseif cnt == 2
    tmp = γ * (2*uprev - 1//2*uprev2 + dt*(2*du₂ - k₁))
  elseif cnt == 3
    tmp = γ * (3*uprev - 3//2*uprev2 + 1//3*uprev3 + dt*(3*(du₂ - k₁) + k₂))
  else
    tmp = γ * (4*uprev - 3*uprev2 + 4//3*uprev3 - 1//4*uprev4 + dt*(4*du₂ - 6*k₁ + 4*k₂ - k₃))
  end
  nlsolver.tmp = tmp

  # Implicit part
  # precalculations
  γdt = γ*dt
  markfirststage!(nlsolver)

  # initial guess
  if alg.extrapolant == :linear
    z = dt*du₁
  else # :constant
    z = zero(u)
  end
  nlsolver.z = z

  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  u = nlsolver.tmp + γ*z

  cnt == 4 && ( cache.uprev4 = uprev3; cache.k₃ = k₂ )
  cnt >= 3 && ( cache.uprev3 = uprev2; cache.k₂ = k₁ )
              ( cache.uprev2 = uprev;  cache.k₁ = du₂ )
  cache.du₁ = f1(u, p, t+dt)
  cache.du₂ = f2(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.destats.nf2 += 1
  integrator.fsallast = cache.du₁ + cache.du₂
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::SBDFCache)
  @unpack uprev, p, t = integrator
  @unpack f1, f2 = integrator.f
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f1(cache.du₁, uprev, p, t)
  f2(cache.du₂, uprev, p, t)
  integrator.destats.nf += 1
  integrator.destats.nf2 += 1
  @.. integrator.fsalfirst = cache.du₁ + cache.du₂
end

function perform_step!(integrator, cache::SBDFCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,alg = integrator
  @unpack uprev2,uprev3,uprev4,k₁,k₂,k₃,du₁,du₂,nlsolver = cache
  @unpack tmp,z = nlsolver
  @unpack f1, f2 = integrator.f
  cnt = cache.cnt = min(alg.order, integrator.iter+1)
  integrator.iter == 1 && !integrator.u_modified && ( cnt = cache.cnt = 1 )
  nlsolver.γ = γ = inv(γₖ[cnt])
  # Explicit part
  if cnt == 1
    @.. tmp = uprev + dt*du₂
  elseif cnt == 2
    @.. tmp = γ * (2*uprev - 1//2*uprev2 + dt*(2*du₂ - k₁))
  elseif cnt == 3
    @.. tmp = γ * (3*uprev - 3//2*uprev2 + 1//3*uprev3 + dt*(3*(du₂ - k₁) + k₂))
  else
    @.. tmp = γ * (4*uprev - 3*uprev2 + 4//3*uprev3 - 1//4*uprev4 + dt*(4*du₂ - 6*k₁ + 4*k₂ - k₃))
  end
  # Implicit part
  # precalculations
  γdt = γ*dt
  markfirststage!(nlsolver)

  # initial guess
  if alg.extrapolant == :linear
    @.. z = dt*du₁
  else # :constant
    @.. z = zero(eltype(u))
  end

  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  @.. u = tmp + γ*z

  cnt == 4 && ( cache.uprev4 .= uprev3; cache.k₃ .= k₂ )
  cnt >= 3 && ( cache.uprev3 .= uprev2; cache.k₂ .= k₁ )
              ( cache.uprev2 .= uprev;  cache.k₁ .= du₂ )
  f1(du₁, u, p, t+dt)
  f2(du₂, u, p, t+dt)
  integrator.destats.nf += 1
  integrator.destats.nf2 += 1
  @.. integrator.fsallast = du₁ + du₂
end

# QNDF1

function initialize!(integrator, cache::QNDF1ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::QNDF1ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uprev2,D,D2,R,U,dtₙ₋₁,nlsolver = cache
  alg = unwrap_alg(integrator, true)
  cnt = integrator.iter
  k = 1
  if cnt == 1
    κ = zero(alg.kappa)
  else
    κ = alg.kappa
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
  nlsolver.tmp = u₀ - ϕ

  γdt = γ*dt
  markfirststage!(nlsolver)

  # initial guess
  nlsolver.z = dt*integrator.fsalfirst
  nlsolver.γ = γ

  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  u = nlsolver.tmp + γ*z

  if integrator.opts.adaptive && integrator.success_iter > 0
    D2[1] = u - uprev
    D2[2] = D2[1] - D[1]
    utilde = (κ*γ₁ + inv(k+1)) * D2[2]
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
    if integrator.EEst > one(integrator.EEst)
      return
    end
  else
    integrator.EEst = one(integrator.EEst)
  end
  cache.dtₙ₋₁ = dt
  cache.uprev2 = uprev
  integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::QNDF1Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

function perform_step!(integrator,cache::QNDF1Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uprev2,D,D2,R,U,dtₙ₋₁,utilde,atmp,nlsolver = cache
  @unpack tmp,z = nlsolver
  alg = unwrap_alg(integrator, true)
  cnt = integrator.iter
  k = 1
  if cnt == 1
    κ = zero(alg.kappa)
  else
    κ = alg.kappa
    ρ = dt/dtₙ₋₁
    @.. D[1] = uprev - uprev2 # backward diff
    if ρ != 1
      R!(k,ρ,cache)
      @.. D[1] = D[1] * (R[1] * U[1])
    end
  end

  # precalculations
  γ₁ = 1//1
  nlsolver.γ = γ = inv((1-κ)*γ₁)
  @.. tmp = uprev + D[1] - γ * (γ₁*D[1])

  γdt = γ*dt
  markfirststage!(nlsolver)

  # initial guess
  @.. z = dt*integrator.fsalfirst

  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  @.. u = tmp + γ*z

  if integrator.opts.adaptive && integrator.success_iter > 0
    @.. D2[1] = u - uprev
    @.. D2[2] = D2[1] - D[1]
    @.. utilde = (κ*γ₁ + inv(k+1)) * D2[2]
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
    if integrator.EEst > one(integrator.EEst)
      return
    end
  else
    integrator.EEst = one(integrator.EEst)
  end
  cache.dtₙ₋₁ = dt
  cache.uprev2 .= uprev
  f(integrator.fsallast, u, p, t+dt)
  integrator.destats.nf += 1
end

function initialize!(integrator, cache::QNDF2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::QNDF2ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uprev2,uprev3,dtₙ₋₁,dtₙ₋₂,D,D2,R,U,nlsolver = cache
  alg = unwrap_alg(integrator, true)
  cnt = integrator.iter
  k = 2
  if cnt == 1 || cnt == 2
    κ = zero(alg.kappa)
    γ₁ = 1//1
    γ₂ = 1//1
  elseif dtₙ₋₁ != dtₙ₋₂
    κ = alg.kappa
    γ₁ = 1//1
    γ₂ = 1//1 + 1//2
    ρ₁ = dt/dtₙ₋₁
    ρ₂ = dt/dtₙ₋₂
    D[1] = uprev - uprev2
    D[1] = D[1] * ρ₁
    D[2] = D[1] - ((uprev2 - uprev3) * ρ₂)
  else
    κ = alg.kappa
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
  nlsolver.γ = γ = inv((1-κ)*γ₂)
  u₀ = uprev + D[1] + D[2]
  ϕ = γ * (γ₁*D[1] + γ₂*D[2])
  nlsolver.tmp = u₀ - ϕ

  γdt = γ*dt
  markfirststage!(nlsolver)

  # initial guess
  nlsolver.z = dt*integrator.fsalfirst

  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  u = nlsolver.tmp + γ*z

  if integrator.opts.adaptive
    if integrator.success_iter == 0
      integrator.EEst = one(integrator.EEst)
    elseif integrator.success_iter == 1
      utilde = (u - uprev) - ((uprev - uprev2) * dt/dtₙ₋₁)
      atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      integrator.EEst = integrator.opts.internalnorm(atmp,t)
    else
      D2[1] = u - uprev
      D2[2] = D2[1] - D[1]
      D2[3] = D2[2] - D[2]
      utilde = (κ*γ₂ + inv(k+1)) * D2[3]
      atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      integrator.EEst = integrator.opts.internalnorm(atmp,t)
    end
  end
  if integrator.EEst > one(integrator.EEst)
    return
  end

  cache.uprev3 = uprev2
  cache.uprev2 = uprev
  cache.dtₙ₋₂ = dtₙ₋₁
  cache.dtₙ₋₁ = dt
  integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
  return
end

function initialize!(integrator, cache::QNDF2Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

function perform_step!(integrator,cache::QNDF2Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uprev2,uprev3,dtₙ₋₁,dtₙ₋₂,D,D2,R,U,utilde,atmp,nlsolver = cache
  tmp = nlsolver.tmp
  alg = unwrap_alg(integrator, true)
  cnt = integrator.iter
  k = 2
  if cnt == 1 || cnt == 2
    κ = zero(alg.kappa)
    γ₁ = 1//1
    γ₂ = 1//1
  elseif dtₙ₋₁ != dtₙ₋₂
    κ = alg.kappa
    γ₁ = 1//1
    γ₂ = 1//1 + 1//2
    ρ₁ = dt/dtₙ₋₁
    ρ₂ = dt/dtₙ₋₂
    @.. D[1] = uprev - uprev2
    @.. D[1] = D[1] * ρ₁
    @.. D[2] = D[1] - ((uprev2 - uprev3) * ρ₂)
  else
    κ = alg.kappa
    γ₁ = 1//1
    γ₂ = 1//1 + 1//2
    ρ = dt/dtₙ₋₁
    # backward diff
    @.. D[1] = uprev - uprev2
    @.. D[2] = D[1] - (uprev2 - uprev3)
    if ρ != 1
      R!(k,ρ,cache)
      R .= R * U
      @.. D[1] = D[1] * R[1,1] + D[2] * R[2,1]
      @.. D[2] = D[1] * R[1,2] + D[2] * R[2,2]
    end
  end

  # precalculations
  nlsolver.γ = γ = inv((1-κ)*γ₂)
  @.. tmp = uprev + D[1] + D[2] - γ * (γ₁*D[1] + γ₂*D[2])

  γdt = γ*dt
  markfirststage!(nlsolver)

  # initial guess
  @.. nlsolver.z = dt*integrator.fsalfirst

  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  @.. u = tmp + γ*z

  if integrator.opts.adaptive
    if integrator.success_iter == 0
      integrator.EEst = one(integrator.EEst)
    elseif integrator.success_iter == 1
      @.. utilde = (u - uprev) - ((uprev - uprev2) * dt/dtₙ₋₁)
      calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      integrator.EEst = integrator.opts.internalnorm(atmp,t)
    else
      @.. D2[1] = u - uprev
      @.. D2[2] = D2[1] - D[1]
      @.. D2[3] = D2[2] - D[2]
      @.. utilde = (κ*γ₂ + inv(k+1)) * D2[3]
      calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      integrator.EEst = integrator.opts.internalnorm(atmp,t)
    end
  end
  if integrator.EEst > one(integrator.EEst)
    return
  end

  cache.uprev3 .= uprev2
  cache.uprev2 .= uprev
  cache.dtₙ₋₂ = dtₙ₋₁
  cache.dtₙ₋₁ = dt
  f(integrator.fsallast, u, p, t+dt)
  integrator.destats.nf += 1
  return
end

function initialize!(integrator, cache::QNDFConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::QNDFConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack udiff,dts,order,max_order,D,D2,R,U,nlsolver = cache
  k = order
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
      n = k+1
      if cnt == 3
        n = k
      end
      for i = 1:n
        D2[1,i] = udiff[i]
      end
      backward_diff!(cache,D,D2,k)
      if ρ != 1
        U!(k,U)
        R!(k,ρ,cache)
        R .= R * U
        reinterpolate_history!(cache,D,R,k)
      end
    else
      n = k+1
      if cnt == 3
        n = k
      end
      for i = 1:n
        D2[1,i] = udiff[i] * dt/dts[i]
      end
      backward_diff!(cache,D,D2,k)
    end
  else
    γ = 1//1
  end
  nlsolver.γ = γ
  # precalculations
  u₀ = uprev + sum(D)  # u₀ is predicted value
  ϕ = zero(γ)
  for i = 1:k
    ϕ += γₖ[i]*D[i]
  end
  ϕ *= γ
  nlsolver.tmp = u₀ - ϕ
  γdt = γ*dt
  markfirststage!(nlsolver)
  # initial guess
  nlsolver.z = dt*integrator.fsalfirst

  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  u = nlsolver.tmp + γ*z

  if integrator.opts.adaptive
    if cnt == 1
      integrator.EEst = one(integrator.EEst)
    elseif cnt == 2
      utilde = (u - uprev) - (udiff[1] * dt/dts[1])
      atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      integrator.EEst = integrator.opts.internalnorm(atmp,t)
    else
      δ = u - uprev
      for i = 1:k
        δ -= D[i]
      end
      utilde = (κ*γₖ[k] + inv(k+1)) * δ
      atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      integrator.EEst = integrator.opts.internalnorm(atmp,t)
    end

    if cnt == 1
      cache.order = 1
    elseif cnt <= 3
      cache.order = 2
    else
      errm1 = 0
      if k > 1
        utildem1 = (κ*γₖ[k-1] + inv(k)) * D[k]
        atmpm1 = calculate_residuals(utildem1, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
        errm1 = integrator.opts.internalnorm(atmpm1, t)
      end
      backward_diff!(cache,D,D2,k+1,false)
      δ = u - uprev
      for i = 1:(k+1)
        δ -= D2[i,1]
      end
      utildep1 = (κ*γₖ[k+1] + inv(k+2)) * δ
      atmpp1 = calculate_residuals(utildep1, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      errp1 = integrator.opts.internalnorm(atmpp1, t)
      pass = stepsize_and_order!(cache, integrator.EEst, errm1, errp1, dt, k)
      if pass == false
        cache.c = cache.c + 1
        fill!(D, zero(u)); fill!(D2, zero(u))
        fill!(R, zero(t)); fill!(U, zero(t))
        return
      end
      cache.c = 0
    end # cnt == 1
  end # integrator.opts.adaptive
  for i = 6:-1:2
    dts[i] = dts[i-1]
    udiff[i] = udiff[i-1]
  end
  dts[1] = dt
  udiff[1] = u - uprev
  fill!(D, zero(u)); fill!(D2, zero(u))
  fill!(R, zero(t)); fill!(U, zero(t))

  integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::QNDFCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

function perform_step!(integrator,cache::QNDFCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack udiff,dts,order,max_order,D,D2,R,U,utilde,atmp,nlsolver = cache
  tmp = nlsolver.tmp
  cnt = integrator.iter
  k = order
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
      n = k+1
      if cnt == 3
        n = k
      end
      for i = 1:n
        D2[1,i] .= udiff[i]
      end
      backward_diff!(cache,D,D2,k)
      if ρ != 1
        U!(k,U)
        R!(k,ρ,cache)
        R .= R * U
        reinterpolate_history!(cache,D,R,k)
      end
    else
      n = k+1
      if cnt == 3
        n = k
      end
      for i = 1:n
        @.. D2[1,i] = udiff[i] * dt/dts[i]
      end
      backward_diff!(cache,D,D2,k)
    end
  else
    γ = one(γ)
  end
  nlsolver.γ = γ
  # precalculations
  ϕ = fill!(utilde, zero(eltype(u)))
  for i = 1:k
    @.. ϕ += γₖ[i]*D[i]
  end
  @.. ϕ *= γ
  tm = fill!(nlsolver.ztmp, zero(eltype(u)))
  for i = 1:k
    @.. tm += D[i]
  end
  @.. nlsolver.tmp = uprev + tm - ϕ

  γdt = γ*dt
  markfirststage!(nlsolver)
  # initial guess
  @.. nlsolver.z = dt*integrator.fsalfirst

  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  @.. u = nlsolver.tmp + γ*z


  if integrator.opts.adaptive
    if cnt == 1
      integrator.EEst = one(integrator.EEst)
    elseif cnt == 2
      @.. utilde = (u - uprev) - (udiff[1] * dt/dts[1])
      calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      integrator.EEst = integrator.opts.internalnorm(atmp,t)
    else
      @.. tmp = u - uprev
      for i = 1:k
        @.. tmp -= D[i]
      end
      @.. utilde = (κ*γₖ[k] + inv(k+1)) * tmp
      calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      integrator.EEst = integrator.opts.internalnorm(atmp,t)
    end

    if cnt == 1
      cache.order = 1
    elseif cnt <= 3
      cache.order = 2
    else
      errm1 = 0
      if k > 1
        @.. utilde = (κ*γₖ[k-1] + inv(k)) * D[k]
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
        errm1 = integrator.opts.internalnorm(atmp,t)
      end
      backward_diff!(cache,D,D2,k+1,false)
      @.. tmp = u - uprev
      for i = 1:(k+1)
        @.. tmp -= D2[i,1]
      end
      @.. utilde = (κ*γₖ[k+1] + inv(k+2)) * tmp
      calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      errp1 = integrator.opts.internalnorm(atmp,t)
      pass = stepsize_and_order!(cache, integrator.EEst, errm1, errp1, dt, k)
      if pass == false
        for i = 1:5
          fill!(D[i], zero(eltype(u)))
        end
        for i = 1:6
          for j = 1:6
            fill!(D2[i,j], zero(eltype(u)))
          end
        end
        fill!(R, zero(t)); fill!(U, zero(t))
        cache.c = cache.c + 1
        return
      end
      cache.c = 0
    end # cnt == 1
  end # integrator.opts.adaptive
  swap_tmp = udiff[6]
  for i = 6:-1:2
    dts[i] = dts[i-1]
    udiff[i] = udiff[i-1]
  end
  udiff[1] = swap_tmp
  dts[1] = dt
  @.. udiff[1] = u - uprev
  for i = 1:5
    fill!(D[i], zero(eltype(u)))
  end
  for i = 1:6, j = 1:6
    fill!(D2[i,j], zero(eltype(u)))
  end
  fill!(R, zero(t)); fill!(U, zero(t))

  f(integrator.fsallast, u, p, t+dt)
  integrator.destats.nf += 1
end

### MEBDF2
function initialize!(integrator, cache::MEBDF2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::MEBDF2ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  nlsolver = cache.nlsolver
  alg = unwrap_alg(integrator, true)
  markfirststage!(nlsolver)

  # initial guess
  if alg.extrapolant == :linear
    nlsolver.z = dt*integrator.fsalfirst
  else # :constant
    nlsolver.z = zero(u)
  end

### STEP 1
  nlsolver.tmp = uprev
  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  z₁ = nlsolver.tmp + z
### STEP 2
  nlsolver.tmp = z₁
  nlsolver.c = 2
  isnewton(nlsolver) && set_new_W!(nlsolver, false)
  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  z₂ = z₁ + z
### STEP 3
  tmp2 = 0.5uprev + z₁ - 0.5z₂
  nlsolver.tmp = tmp2
  nlsolver.c = 1
  isnewton(nlsolver) && set_new_W!(nlsolver, false)
  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  u = tmp2 + z

### finalize
  integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::MEBDF2Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::MEBDF2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z₁,z₂,tmp2,nlsolver = cache
  z = nlsolver.z
  mass_matrix = integrator.f.mass_matrix
  alg = unwrap_alg(integrator, true)
  markfirststage!(nlsolver)

  # initial guess
  if alg.extrapolant == :linear
    @.. z = dt*integrator.fsalfirst
  else # :constant
    z .= zero(eltype(u))
  end

### STEP 1
 nlsolver.tmp = uprev
 z = nlsolve!(nlsolver, integrator, cache, repeat_step)
 nlsolvefail(nlsolver) && return
 @.. z₁ = uprev + z
### STEP 2
 nlsolver.tmp = z₁
 nlsolver.c = 2
 isnewton(nlsolver) && set_new_W!(nlsolver, false)
 z = nlsolve!(nlsolver, integrator, cache, repeat_step)
 nlsolvefail(nlsolver) && return
 @.. z₂ = z₁ + z
### STEP 3
 # z .= zero(eltype(u))
 @.. tmp2 = 0.5uprev + z₁ - 0.5z₂
 nlsolver.tmp = tmp2
 nlsolver.c = 1
 isnewton(nlsolver) && set_new_W!(nlsolver, false)
 z = nlsolve!(nlsolver, integrator, cache, repeat_step)
 nlsolvefail(nlsolver) && return
 @.. u = tmp2 + z

### finalize
 f(integrator.fsallast,u,p,t+dt)
 integrator.destats.nf += 1
end
