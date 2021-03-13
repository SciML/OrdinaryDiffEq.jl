function initialize!(integrator, cache::Union{DImplicitEulerCache})
  integrator.kshortsize = 2
  @unpack k₁,k₂ = cache
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [k₁,k₂]
  integrator.k[1] .= integrator.du
  nothing
end

@muladd function perform_step!(integrator, cache::DImplicitEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  alg = unwrap_alg(integrator, true)
  @unpack nlsolver = cache

  nlsolver.z = zero(u)
  nlsolver.tmp = zero(u)
  nlsolver.γ = 1
  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  u = uprev + z

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

    tmp = r*integrator.opts.internalnorm.((u - uprev)/dt1 - (uprev - uprev2)/dt2,t)
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  else
    integrator.EEst = 1
  end

  integrator.u = u
  integrator.du = (u-uprev)/dt

  if integrator.opts.calck
    integrator.k[2] = integrator.k[1]
    integrator.k[1] = integrator.du
  end
end


@muladd function perform_step!(integrator, cache::DImplicitEulerCache, repeat_step=false)
  @unpack t,dt,uprev,du,u,f,p = integrator
  @unpack atmp,nlsolver = cache
  @unpack tmp = nlsolver
  alg = unwrap_alg(integrator, true)

  @. nlsolver.z = false
  @. nlsolver.tmp = false
  nlsolver.γ = 1
  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  @.. u = uprev + z
  @.. du = z * inv(dt)

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

    @.. tmp = r*integrator.opts.internalnorm((u - uprev)/dt1 - (uprev - uprev2)/dt2,t)
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  else
    integrator.EEst = 1
  end

  if integrator.opts.calck
    integrator.k[2] .= integrator.k[1]
    integrator.k[1] .= du
  end
end

function initialize!(integrator, cache::DABDF2ConstantCache) end

@muladd function perform_step!(integrator, cache::DABDF2ConstantCache, repeat_step=false)
  @unpack t,f,p = integrator
  @unpack dtₙ₋₁,nlsolver = cache
  alg = unwrap_alg(integrator, true)
  dtₙ, uₙ, uₙ₋₁, uₙ₋₂ = integrator.dt, integrator.u, integrator.uprev, integrator.uprev2

  if integrator.iter == 1 && !integrator.u_modified
    cache.dtₙ₋₁ = dtₙ
    perform_step!(integrator, cache.eulercache, repeat_step)
    integrator.fsalfirst = @.. (integrator.u - integrator.uprev) / dtₙ
    cache.fsalfirstprev = integrator.fsalfirst
    return
  end

  # precalculations
  ρ = dtₙ/dtₙ₋₁
  c1 = ρ^2/(1+2ρ)

  nlsolver.γ = (1+ρ)/(1+2ρ)
  nlsolver.α = 1//1

  nlsolver.z = zero(uₙ)

  nlsolver.tmp = -c1 * uₙ₋₁ + c1 * uₙ₋₂
  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return

  uₙ = uₙ₋₁ + z
  integrator.fsallast = @.. z/dtₙ

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

  integrator.u = uₙ
  integrator.du = du = (nlsolver.α * z + nlsolver.tmp) * inv(nlsolver.γ * dtₙ)

  if integrator.opts.calck
    integrator.k[2] = integrator.k[1]
    integrator.k[1] = integrator.du
  end
  return
end

function initialize!(integrator, cache::DABDF2Cache)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)

  integrator.kshortsize = 2
  @unpack k₁,k₂ = cache.eulercache
  resize!(integrator.k, integrator.kshortsize)
  integrator.k .= [k₁,k₂]
  integrator.k[1] .= integrator.du
  nothing
end

@muladd function perform_step!(integrator, cache::DABDF2Cache, repeat_step=false)
  @unpack t,dt,du,f,p = integrator
  @unpack atmp,dtₙ₋₁,nlsolver = cache
  @unpack z,tmp = nlsolver
  alg = unwrap_alg(integrator, true)
  uₙ,uₙ₋₁,uₙ₋₂,dtₙ = integrator.u,integrator.uprev,integrator.uprev2,integrator.dt

  if integrator.iter == 1 && !integrator.u_modified
    cache.dtₙ₋₁ = dtₙ
    perform_step!(integrator, cache.eulercache, repeat_step)
    @.. integrator.fsalfirst = (uₙ - uₙ₋₁) / dt
    cache.fsalfirstprev .= integrator.fsalfirst
    return
  end

  # precalculations
  ρ = dtₙ/dtₙ₋₁
  c1 = ρ^2/(1+2ρ)

  nlsolver.γ = (1+ρ)/(1+2ρ)
  nlsolver.α = 1//1
  @.. nlsolver.tmp = -c1 * uₙ₋₁ + c1 * uₙ₋₂
  nlsolver.z .= zero(eltype(z))
  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return

  @.. uₙ = uₙ₋₁ + z
  @.. du = (nlsolver.α * z + nlsolver.tmp) * inv(nlsolver.γ * dt)

  @.. integrator.fsallast = du
  integrator.destats.nf += 1
  if integrator.opts.adaptive
    btilde0 = (dtₙ₋₁+dtₙ)*1//6
    btilde1 = 1+ρ
    btilde2 = ρ
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
