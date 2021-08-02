function initialize!(integrator, cache::DImplicitEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.k[1] = integrator.du
end

function initialize!(integrator, cache::DImplicitEulerCache)
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

function initialize!(integrator, cache::DABDF2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.k[1] = integrator.du
end

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

function initialize!(integrator, cache::DFBDFConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.k[1] = integrator.du
end

function perform_step!(integrator, cache::DFBDFConstantCache{max_order}, repeat_step=false) where max_order
  @unpack ts,u_history,order,u_corrector,bdf_coeffs,r,nlsolver,weights,ts_tmp = cache
  @unpack t,dt,u,f,p,uprev = integrator

  if integrator.u_modified
    order = 1
    cache.consfailcnt = cache.nconsteps = 0
    fill!(weights,zero(eltype(weights)))
    fill!(ts,zero(eltype(ts)))
    fill!(u_history,zero(eltype(u_history)))
    fill!(u_corrector,zero(eltype(u_corrector)))
    cache.nonevesuccsteps = 0
  end
  @unpack nonevesuccsteps,consfailcnt,nconsteps = cache

  k = order
  if nonevesuccsteps == 0
    weights[1] = 1/dt
    ts[1] = t
    @.. u_history[:,1] = $_vec(uprev)
  elseif nonevesuccsteps == 1
    weights[1] = inv(t-ts[1])
    weights[2] = inv(ts[1]-t)
    ts[2] = ts[1]
    ts[1] = t
    @.. @views u_history[:,2] = u_history[:,1]
    @.. u_history[:,1] = $_vec(uprev)
  elseif consfailcnt == 0
    for i in k+2:-1:2
      ts[i] = ts[i-1]
      @.. @views u_history[:,i] = u_history[:,i-1]
    end
    ts[1] = t
    @.. u_history[:,1] = $_vec(uprev)
  end
  
  if nonevesuccsteps >= 1
    compute_weights!(ts,k,weights)
  end
    
  u₀ = zero(u)
  if nonevesuccsteps >= 1
    u₀ = calc_Lagrange_interp(k,weights,t+dt,ts,u_history,u₀)
  else
    u₀ = u
  end
  markfirststage!(nlsolver)
  
  nlsolver.z = zero(u₀)

  equi_ts = zeros(k-1)
  for i in 1:k-1
    equi_ts[i] = t - dt*i
  end

  fill!(u_corrector,zero(eltype(u)))
  if u isa Number
    for i in 1:k-1
      u_corrector[i] = calc_Lagrange_interp(k,weights,equi_ts[i],ts,u_history,u_corrector[i])
    end
    tmp = uprev * bdf_coeffs[k,2]
    for i in 1:k-1
      tmp += u_corrector[i] * bdf_coeffs[k,i+2]
    end
  else
    for i in 1:k-1
      @.. @views u_corrector[:,i] = $calc_Lagrange_interp(k,weights,equi_ts[i],ts,u_history,u_corrector[:,i])
    end
    tmp = uprev * bdf_coeffs[k,2]
    vc = _vec(tmp)
    for i in 1:k-1
      @.. @views vc += u_corrector[:,i] * bdf_coeffs[k,i+2]
    end
  end

  nlsolver.tmp = tmp + u₀
  β₀ = bdf_coeffs[k,1]
  α₀ = 1//1 #bdf_coeffs[k,1]
  nlsolver.γ = β₀
  nlsolver.α = α₀
  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) &&  return

  u = u₀ + z

  for j in 2:k
    r[j] = (1-j)
    for i in 2:k+1
      r[j] *= ((t+dt-j*dt)-ts[i])/(i*dt) #TODO: This should be noticed that whether it uses the correct ts elements.
    end
  end

  terkp1 = z
  for j in 1:k+1
    terkp1 *= j*dt/(t+dt-ts[j])
  end

  lte = -1/(1+k)
  for j in 2:k
    lte -= bdf_coeffs[k,j]*r[j]
  end
  lte *= terkp1

  if integrator.opts.adaptive
    for i in 1:k+1
      ts_tmp[i+1] = ts[i]
    end
    ts_tmp[1] = t+dt
    atmp = calculate_residuals(_vec(lte), _vec(uprev), _vec(u), integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)

    fd_weights = calc_finite_difference_weights(ts_tmp,t+dt,k,Val(max_order))
    terk = @.. fd_weights[1,k+1] * u

    if u isa Number
      for i in 2:k+1
        terk += fd_weights[i,k+1] * u_history[i-1]
      end
      terk *= abs(dt^(k))
    else
      vc = _vec(terk)
      for i in 2:k+1
        @.. @views vc += fd_weights[i,k+1] * u_history[:,i-1]
      end
      terk *= abs(dt^(k))
    end

    atmp = calculate_residuals(_vec(terk), _vec(uprev), _vec(u), integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    cache.terk = integrator.opts.internalnorm(atmp,t)
    
    if k > 1
      fd_weights = calc_finite_difference_weights(ts_tmp,t+dt,k-1,Val(max_order))
      terkm1 = fd_weights[1,k] * u

      if u isa Number
        for i in 2:k
          terkm1 += fd_weights[i,k] * u_history[i-1]
        end
        terkm1 *= abs(dt^(k-1))
      else
        vc = _vec(terkm1)
        for i in 2:k
          @.. @views vc += fd_weights[i,k] * u_history[:,i-1]
        end
        terkm1 *= abs(dt^(k-1))
      end
      atmp = calculate_residuals(_vec(terkm1), _vec(uprev), _vec(u), integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      cache.terkm1 = integrator.opts.internalnorm(atmp,t)
    end
    if k > 2
      fd_weights = calc_finite_difference_weights(ts_tmp,t+dt,k-2,Val(max_order))
      terkm2 = fd_weights[1,k-1] * u

      if u isa Number
        for i in 2:k-1
          terkm2 += fd_weights[i,k-1] * u_history[i-1]
        end
        terkm2 *= abs(dt^(k-2))
      else
        vc = _vec(terkm2)
        for i in 2:k-1
          @.. @views vc += fd_weights[i,k-1] * u_history[:,i-1]
        end
        terkm2 *= abs(dt^(k-2))
      end
      atmp = calculate_residuals(_vec(terkm2), _vec(uprev), _vec(u), integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      cache.terkm2 = integrator.opts.internalnorm(atmp,t)
    end
    if nconsteps > k+1 && k < max_order
      atmp = calculate_residuals(_vec(terkp1), _vec(uprev), _vec(u), integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      cache.terkp1 = integrator.opts.internalnorm(atmp,t)
    else
      cache.terkp1 = zero(cache.terk)
    end
  end

  integrator.destats.nf += 1
  integrator.u = u
  integrator.fsallast = integrator.du = (nlsolver.α * z + nlsolver.tmp) * inv(nlsolver.γ * dt)
  if integrator.opts.calck
    integrator.k[2] = integrator.k[1]
    integrator.k[1] = integrator.du
  end
end

function initialize!(integrator, cache::DFBDFCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::DFBDFCache{max_order}, repeat_step=false) where max_order
  @unpack ts,u_history,order,u_corrector,bdf_coeffs,r,nlsolver,weights,terk_tmp,terkp1_tmp,atmp,tmp,equi_ts,u₀,ts_tmp = cache
  @unpack t,dt,u,f,p,uprev = integrator

  if integrator.u_modified
    order = 1
    cache.consfailcnt = cache.nconsteps = 0
    fill!(weights,zero(eltype(weights)))
    fill!(ts,zero(eltype(ts)))
    fill!(u_history,zero(eltype(u_history)))
    fill!(u_corrector,zero(eltype(u_corrector)))
    cache.nonevesuccsteps = 0
  end
  @unpack nonevesuccsteps,consfailcnt,nconsteps = cache
  
  k = order
  if nonevesuccsteps == 0
    weights[1] = 1/dt
    ts[1] = t
    @.. u_history[:,1] = $_vec(uprev)
  elseif nonevesuccsteps == 1
    weights[1] = inv(t-ts[1])
    weights[2] = inv(ts[1]-t)
    ts[2] = ts[1]
    ts[1] = t
    @.. @views u_history[:,2] = u_history[:,1]
    @.. u_history[:,1] = $_vec(uprev)
  elseif consfailcnt == 0 && nlsolver.status == Convergence
    for i in k+2:-1:2
      ts[i] = ts[i-1]
      @.. @views u_history[:,i] = u_history[:,i-1]
    end
    ts[1] = t
    @.. u_history[:,1] = $_vec(uprev)
  end
  if nonevesuccsteps >= 1
    compute_weights!(ts,k,weights)
  end
    
  @.. u₀ = zero(u)
  if nonevesuccsteps >= 1
    calc_Lagrange_interp!(k,weights,t+dt,ts,u_history,u₀)
  else
    @.. u₀ = u
  end
  markfirststage!(nlsolver)
  
  for i in 1:k-1
    equi_ts[i] = t - dt*i
  end

  fill!(u_corrector,zero(eltype(u)))
  for i in 1:k-1
    @views calc_Lagrange_interp!(k,weights,equi_ts[i],ts,u_history,u_corrector[:,i])
  end

  @.. tmp = uprev * bdf_coeffs[k,2]
  vc = _vec(tmp)
  for i in 1:k-1
    @.. @views vc += u_corrector[:,i] * bdf_coeffs[k,i+2]
  end

  @.. nlsolver.tmp = tmp + u₀
  @.. nlsolver.z = zero(eltype(nlsolver.z))#newton dae uprev!=u₀
  nlsolver.γ = bdf_coeffs[k,1]
  nlsolver.α = 1//1
  z = nlsolve!(nlsolver, integrator, cache, repeat_step)
  nlsolvefail(nlsolver) && return
  @.. u = z + u₀

  for j in 2:k
    r[j] = (1-j)
    for i in 2:k+1
      r[j] *= ((t+dt-j*dt)-ts[i])/(i*dt) #TODO: This should be noticed that whether it uses the correct ts elements.
    end
  end

  @.. terkp1_tmp = z
  for j in 1:k+1
    @.. terkp1_tmp *= j*dt/(t+dt-ts[j])
  end

  lte = -1/(1+k)
  for j in 2:k
    lte -= bdf_coeffs[k,j]*r[j]
  end
  @.. terk_tmp = lte * terkp1_tmp
  if integrator.opts.adaptive
    @unpack abstol, reltol, internalnorm = integrator.opts
    for i in 1:k+1
      ts_tmp[i+1] = ts[i]
    end
    ts_tmp[1] = t+dt
    calculate_residuals!(atmp, _vec(terk_tmp), _vec(uprev), _vec(u), abstol, reltol, internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
    fd_weights = calc_finite_difference_weights(ts_tmp,t+dt,k,Val(max_order))
    @.. terk_tmp = fd_weights[1,k+1] * u
    vc = _vec(terk_tmp)
    for i in 2:k+1
      @.. @views vc += fd_weights[i,k+1] * u_history[:,i-1]
    end
    @.. terk_tmp *= abs(dt^(k))
    calculate_residuals!(atmp, _vec(terk_tmp), _vec(uprev), _vec(u), abstol, reltol, internalnorm, t)
    cache.terk = integrator.opts.internalnorm(atmp,t)
    
    if k > 1
      fd_weights = calc_finite_difference_weights(ts_tmp,t+dt,k-1,Val(max_order))
      @.. terk_tmp = fd_weights[1,k] * u
      vc = _vec(terk_tmp)
      for i in 2:k
        @.. @views vc += fd_weights[i,k] *u_history[:,i-1]
      end
      @.. terk_tmp *= abs(dt^(k-1))
      calculate_residuals!(atmp, _vec(terk_tmp), _vec(uprev), _vec(u), integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      cache.terkm1 = integrator.opts.internalnorm(atmp,t)
    end
    if k > 2
      fd_weights = calc_finite_difference_weights(ts_tmp,t+dt,k-2,Val(max_order))
      @.. terk_tmp = fd_weights[1,k-1] * u
      vc = _vec(terk_tmp)
      for i in 2:k-1
        @.. @views vc += fd_weights[i,k-1] *u_history[:,i-1]
      end
      @.. terk_tmp *= abs(dt^(k-2))
      calculate_residuals!(atmp, _vec(terk_tmp), _vec(uprev), _vec(u), integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      cache.terkm2 = integrator.opts.internalnorm(atmp,t)
    end
    if nconsteps > k+1 && k < max_order
      calculate_residuals!(atmp, _vec(terkp1_tmp), _vec(uprev), _vec(u), integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      cache.terkp1 = integrator.opts.internalnorm(atmp,t)
    else
      cache.terkp1 = zero(cache.terkp1)
    end
  end
  integrator.destats.nf += 1
  @.. integrator.fsallast = integrator.du = (nlsolver.α * z + nlsolver.tmp) * inv(nlsolver.γ * dt)
end
