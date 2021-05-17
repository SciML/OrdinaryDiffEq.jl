function initialize!(integrator, cache::ROCK2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  cache.max_stage = (integrator.alg.max_stages < 1 || integrator.alg.max_stages > 200) ? 200 : integrator.alg.max_stages
  cache.min_stage = (integrator.alg.min_stages > cache.max_stage) ? cache.max_stage : integrator.alg.min_stages
  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ROCK2ConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack ms, fp1, fp2, recf = cache
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
  # The the number of degree for Chebyshev polynomial
  mdeg = Int(floor(sqrt((1.5 + abs(dt)*integrator.eigen_est)/0.811) + 1))
  mdeg = min(max(mdeg,cache.min_stage), cache.max_stage)
  cache.mdeg = max(mdeg, 3) - 2
  choosedeg!(cache)
  # recurrence
  # for the first stage
  tᵢ₋₁ = t + dt*recf[cache.start]
  tᵢ₋₂ = t + dt*recf[cache.start]
  tᵢ₋₃ = t
  uᵢ₋₂ = copy(uprev)
  uᵢ₋₁ = uprev + (dt*recf[cache.start])*fsalfirst
  cache.mdeg < 2 && ( u = uᵢ₋₁ )
  # for the second to the ms[cache.mdeg] th stages
  for i in 2:cache.mdeg
    μ, κ = recf[cache.start + (i - 2)*2 + 1], recf[cache.start + (i - 2)*2 + 2]
    ν = -1 - κ
    u = f(uᵢ₋₁, p, tᵢ₋₁)
    tᵢ₋₁ = dt*μ - ν*tᵢ₋₂ - κ*tᵢ₋₃
    u = (dt*μ)*u - ν*uᵢ₋₁ - κ*uᵢ₋₂
    i < cache.mdeg && (uᵢ₋₂ = uᵢ₋₁; uᵢ₋₁ = u)
    tᵢ₋₃ = tᵢ₋₂
    tᵢ₋₂ = tᵢ₋₁
  end # end if
  # two-stage finishing procedure.
  δt₁ = dt*fp1[cache.deg_index]
  δt₂ = dt*fp2[cache.deg_index]
  uᵢ₋₂  = f(u, p, tᵢ₋₁)
  integrator.destats.nf += 1
  uᵢ₋₁ = u + δt₁*uᵢ₋₂
  tᵢ₋₁ += δt₁
  u = f(uᵢ₋₁, p, tᵢ₋₁)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    tmp = δt₂*(u - uᵢ₋₂)
    u = uᵢ₋₁ + δt₁*u + tmp
  else
    u = uᵢ₋₁ + δt₁*u + δt₂*(u - uᵢ₋₂)
  end
  # error estimate
  if integrator.opts.adaptive
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.u = u
end

function initialize!(integrator, cache::ROCK2Cache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
  integrator.fsallast = cache.k
  cache.constantcache.max_stage = (integrator.alg.max_stages < 1 || integrator.alg.max_stages > 200) ? 200 : integrator.alg.max_stages
  cache.constantcache.min_stage = (integrator.alg.min_stages > cache.constantcache.max_stage) ? cache.constantcache.max_stage : integrator.alg.min_stages

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::ROCK2Cache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack k, tmp, uᵢ₋₂, uᵢ₋₁, atmp = cache
  @unpack ms, fp1, fp2, recf = cache.constantcache
  ccache = cache.constantcache
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
  # The the number of degree for Chebyshev polynomial
  mdeg = Int(floor(sqrt((1.5 + abs(dt)*integrator.eigen_est)/0.811) + 1))
  mdeg = min(max(mdeg,ccache.min_stage), ccache.max_stage)
  ccache.mdeg = max(mdeg, 3) - 2
  choosedeg!(cache)
  # recurrence
  # for the first stage
  tᵢ₋₁ = t + dt*recf[ccache.start]
  tᵢ₋₂ = t + dt*recf[ccache.start]
  tᵢ₋₃ = t
  @.. uᵢ₋₂ = uprev
  @.. uᵢ₋₁ = uprev + (dt*recf[ccache.start])*fsalfirst
  ccache.mdeg < 2 && ( @.. u = uᵢ₋₁ )
  # for the second to the ms[ccache.mdeg] th stages
  for i in 2:ccache.mdeg
    μ, κ = recf[ccache.start + (i - 2)*2 + 1], recf[ccache.start + (i - 2)*2 + 2]
    ν = -1 - κ
    f(k, uᵢ₋₁, p, tᵢ₋₁)
    tᵢ₋₁ = dt*μ - ν*tᵢ₋₂ - κ*tᵢ₋₃
    @.. u = (dt*μ)*k - ν*uᵢ₋₁ - κ*uᵢ₋₂
    if i < ccache.mdeg 
      @.. uᵢ₋₂ = uᵢ₋₁
      @.. uᵢ₋₁ = u
    end
    tᵢ₋₃ = tᵢ₋₂
    tᵢ₋₂ = tᵢ₋₁
  end # end if
  # two-stage finishing procedure.
  δt₁ = dt*fp1[ccache.deg_index]
  δt₂ = dt*fp2[ccache.deg_index]
  f(k, u, p, tᵢ₋₁)
  integrator.destats.nf += 1
  @.. uᵢ₋₁ = u + δt₁*k
  if integrator.opts.adaptive
    @.. tmp = - δt₂*k
  else
    @.. u = -δt₂*k
  end
  c = DiffEqBase.value(sign(δt₁))*integrator.opts.internalnorm(δt₁,t)
  tᵢ₋₁ += c
  f(k, uᵢ₋₁, p, tᵢ₋₁)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    @.. tmp += δt₂*k
    @.. u = uᵢ₋₁ + δt₁*k + tmp
  else
    @.. u += uᵢ₋₁ + (δt₁ + δt₂)*k
  end

  # error estimate
  if integrator.opts.adaptive
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  f(integrator.fsallast, u, p, t+dt)
  integrator.destats.nf += 1
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::ROCK4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  cache.max_stage = (integrator.alg.max_stages < 1 || integrator.alg.max_stages > 152) ? 152 : integrator.alg.max_stages
  cache.min_stage = (integrator.alg.min_stages > cache.max_stage) ? cache.max_stage : integrator.alg.min_stages
  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ROCK4ConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack ms, fpa, fpb, fpbe, recf = cache
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
  # The the number of degree for Chebyshev polynomial
  mdeg = Int(floor(sqrt((3 + abs(dt)*integrator.eigen_est)/0.353) + 1))
  mdeg = min(max(mdeg,cache.min_stage), cache.max_stage)
  cache.mdeg = max(mdeg, 5) - 4
  choosedeg!(cache)
  # recurrence
  # for the first stage
  tᵢ₋₁ = t + dt*recf[cache.start]
  tᵢ₋₂ = t + dt*recf[cache.start]
  tᵢ₋₃ = t
  uᵢ₋₂ = copy(uprev)
  uᵢ₋₁ = uprev + (dt*recf[cache.start])*fsalfirst
  cache.mdeg < 2 && ( u = uᵢ₋₁ )
  # for the second to the cache.mdeg th stages
  for i in 2:cache.mdeg
    μ, κ = recf[cache.start + (i - 2)*2 + 1], recf[cache.start + (i - 2)*2 + 2]
    ν = - 1 - κ
    u = f(uᵢ₋₁, p,tᵢ₋₁)
    tᵢ₋₁ = dt*μ - ν*tᵢ₋₂ - κ*tᵢ₋₃
    u = (dt*μ)*u - ν*uᵢ₋₁ - κ*uᵢ₋₂
    i < cache.mdeg && (uᵢ₋₂ = uᵢ₋₁; uᵢ₋₁ = u)
    tᵢ₋₃ = tᵢ₋₂
    tᵢ₋₂ = tᵢ₋₁
  end

  # These constants correspond to the Buther Tableau coefficients of explicit RK methods
  a₂₁ = dt*fpa[cache.deg_index][1]
  a₃₁ = dt*fpa[cache.deg_index][2]; a₃₂ = dt*fpa[cache.deg_index][3]
  a₄₁ = dt*fpa[cache.deg_index][4]; a₄₂ = dt*fpa[cache.deg_index][5]; a₄₃ = dt*fpa[cache.deg_index][6]
  B₁  = dt*fpb[cache.deg_index][1]; B₂  = dt*fpb[cache.deg_index][2]; B₃  = dt*fpb[cache.deg_index][3]; B₄  = dt*fpb[cache.deg_index][4]
  # coefficients of embedded method for error estimation
  B̂₁  = dt*(fpbe[cache.deg_index][1] - fpb[cache.deg_index][1])
  B̂₂  = dt*(fpbe[cache.deg_index][2] - fpb[cache.deg_index][2])
  B̂₃  = dt*(fpbe[cache.deg_index][3] - fpb[cache.deg_index][3])
  B̂₄  = dt*(fpbe[cache.deg_index][4] - fpb[cache.deg_index][4])
  B̂₅  = dt*fpbe[cache.deg_index][5]

  # 4-stage finishing procedure.
  # Stage-1
  uᵢ₋₁ = f(u, p, tᵢ₋₁)
  integrator.destats.nf += 1
  uᵢ₋₂ = u + a₃₁*uᵢ₋₁
  uᵢ₋₃ = u + a₄₁*uᵢ₋₁
  u    += B₁*uᵢ₋₁
  integrator.opts.adaptive && (tmp = B̂₁*uᵢ₋₁)
  uᵢ₋₁ = u + (a₂₁ - B₁)*uᵢ₋₁

  # Stage-2
  c₂ = a₂₁
  _c₂ = DiffEqBase.value(sign(c₂))*integrator.opts.internalnorm(c₂,t)
  tᵢ₋₂ = tᵢ₋₁ + _c₂
  uᵢ₋₁ = f(uᵢ₋₁, p, tᵢ₋₂)
  integrator.destats.nf += 1
  uᵢ₋₂ += a₃₂*uᵢ₋₁
  uᵢ₋₃ += a₄₂*uᵢ₋₁
  u    += B₂*uᵢ₋₁
  integrator.opts.adaptive && (tmp += B̂₂*uᵢ₋₁)

  # Stage-3
  c₃ = a₃₁ + a₃₂
  _c₃ = DiffEqBase.value(sign(c₃))*integrator.opts.internalnorm(c₃,t)
  tᵢ₋₂ = tᵢ₋₁ + _c₃
  uᵢ₋₂ = f(uᵢ₋₂, p, tᵢ₋₂)
  integrator.destats.nf += 1
  uᵢ₋₃ += a₄₃*uᵢ₋₂
  u    += B₃*uᵢ₋₂
  integrator.opts.adaptive && (tmp += B̂₃*uᵢ₋₂)

  #Stage-4
  c₄ = a₄₁ + a₄₂ + a₄₃
  _c₄ = DiffEqBase.value(sign(c₄))*integrator.opts.internalnorm(c₄,t)
  tᵢ₋₂ = tᵢ₋₁ + _c₄
  uᵢ₋₃ = f(uᵢ₋₃, p, tᵢ₋₂)
  integrator.destats.nf += 1
  u    += B₄*uᵢ₋₃
  integrator.opts.adaptive && (tmp += B̂₄*uᵢ₋₃)

  uᵢ₋₁ = f(u, p, t + dt)
  integrator.destats.nf += 1

  #Error estimate (embedded method of order 3)
  if integrator.opts.adaptive
    tmp  += B̂₅*uᵢ₋₁
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast = uᵢ₋₁
  integrator.u = u
end

function initialize!(integrator, cache::ROCK4Cache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  cache.constantcache.max_stage = (integrator.alg.max_stages < 1 || integrator.alg.max_stages > 152) ? 152 : integrator.alg.max_stages
  cache.constantcache.min_stage = (integrator.alg.min_stages > cache.constantcache.max_stage) ? cache.constantcache.max_stage : integrator.alg.min_stages

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::ROCK4Cache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack uᵢ₋₁, uᵢ₋₂, uᵢ₋₃, tmp, atmp, k= cache
  @unpack ms, fpa, fpb, fpbe, recf = cache.constantcache
  ccache = cache.constantcache
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
  # The the number of degree for Chebyshev polynomial
  mdeg = Int(floor(sqrt((3 + abs(dt)*integrator.eigen_est)/0.353) + 1))
  mdeg = min(max(mdeg,ccache.min_stage), ccache.max_stage)
  ccache.mdeg = max(mdeg, 5) - 4
  choosedeg!(cache)
  # recurrence
  # for the first stage
  tᵢ₋₁ = t + dt*recf[ccache.start]
  tᵢ₋₂ = t + dt*recf[ccache.start]
  tᵢ₋₃ = t
  @.. uᵢ₋₂ = uprev
  @.. uᵢ₋₁ = uprev + (dt*recf[ccache.start])*fsalfirst
  ccache.mdeg < 2 && ( @.. u = uᵢ₋₁ )
  # for the second to the ccache.mdeg th stages
  for i in 2:ccache.mdeg
    μ, κ = recf[ccache.start + (i - 2)*2 + 1], recf[ccache.start + (i - 2)*2 + 2]
    ν = -1 - κ
    f(k, uᵢ₋₁, p, tᵢ₋₁)
    tᵢ₋₁ = (dt*μ) - ν*tᵢ₋₂ - κ*tᵢ₋₃
    @.. u = (dt*μ)*k - ν*uᵢ₋₁ - κ*uᵢ₋₂
    if i < ccache.mdeg 
      @.. uᵢ₋₂ = uᵢ₋₁
      @.. uᵢ₋₁ = u
    end    
    tᵢ₋₃ = tᵢ₋₂
    tᵢ₋₂ = tᵢ₋₁
  end

  # These constants correspond to the Buther Tableau coefficients of explicit RK methods
  a₂₁ = dt*fpa[ccache.deg_index][1]
  a₃₁ = dt*fpa[ccache.deg_index][2]; a₃₂ = dt*fpa[ccache.deg_index][3]
  a₄₁ = dt*fpa[ccache.deg_index][4]; a₄₂ = dt*fpa[ccache.deg_index][5]; a₄₃ = dt*fpa[ccache.deg_index][6]
  B₁  = dt*fpb[ccache.deg_index][1]; B₂  = dt*fpb[ccache.deg_index][2]; B₃  = dt*fpb[ccache.deg_index][3]; B₄  = dt*fpb[ccache.deg_index][4]
  # coefficients of embedded method for error estimation
  B̂₁  = dt*(fpbe[ccache.deg_index][1] - fpb[ccache.deg_index][1])
  B̂₂  = dt*(fpbe[ccache.deg_index][2] - fpb[ccache.deg_index][2])
  B̂₃  = dt*(fpbe[ccache.deg_index][3] - fpb[ccache.deg_index][3])
  B̂₄  = dt*(fpbe[ccache.deg_index][4] - fpb[ccache.deg_index][4])
  B̂₅  = dt*fpbe[ccache.deg_index][5]

  # 4-stage finishing procedure.
  # Stage-1

  f(k, u, p, tᵢ₋₁)
  integrator.destats.nf += 1
  @.. uᵢ₋₂ = u + a₃₁*k
  @.. uᵢ₋₃ = u + a₄₁*k
  @.. uᵢ₋₁ = u + a₂₁*k
  @.. u    += B₁*k
  integrator.opts.adaptive && (@.. tmp = B̂₁*k)

  # Stage-2
  c₂ = a₂₁
  _c₂ = DiffEqBase.value(sign(c₂))*integrator.opts.internalnorm(c₂,t)
  tᵢ₋₂ = tᵢ₋₁ + _c₂
  f(k, uᵢ₋₁, p, tᵢ₋₂)
  integrator.destats.nf += 1
  @.. uᵢ₋₂ += a₃₂*k
  @.. uᵢ₋₃ += a₄₂*k
  @.. u    += B₂*k
  integrator.opts.adaptive && (@.. tmp += B̂₂*k)

  # Stage-3
  c₃ = a₃₁ + a₃₂
  _c₃ = DiffEqBase.value(sign(c₃))*integrator.opts.internalnorm(c₃,t)
  tᵢ₋₂ = tᵢ₋₁ + _c₃
  f(k, uᵢ₋₂, p, tᵢ₋₂)
  integrator.destats.nf += 1
  @.. uᵢ₋₃ += a₄₃*k
  @.. u    += B₃*k
  integrator.opts.adaptive && (@.. tmp += B̂₃*k)

  #Stage-4
  c₄ = a₄₁ + a₄₂ + a₄₃
  _c₄ = DiffEqBase.value(sign(c₄))*integrator.opts.internalnorm(c₄,t)
  tᵢ₋₂ = tᵢ₋₁ + _c₄
  f(k, uᵢ₋₃, p, tᵢ₋₂)
  integrator.destats.nf += 1
  @.. u    += B₄*k
  integrator.opts.adaptive && (tmp += B̂₄*k)

  f(k, u, p, t + dt)
  integrator.destats.nf += 1

  #Error estimate (embedded method of order 3)
  if integrator.opts.adaptive
    tmp  += B̂₅*k
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  @.. integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.destats.nf += 1
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::RKCConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::RKCConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
  # The the number of degree for Chebyshev polynomial
  #maxm = max(2,Int(floor(sqrt(integrator.opts.internalnorm(integrator.opts.reltol,t)/(10*eps(integrator.opts.internalnorm(uprev,t)))))))
  maxm = 50
  mdeg = 1 + Int(floor(sqrt(1.54*abs(dt)*integrator.eigen_est + 1)))
  mdeg = (mdeg > maxm) ? maxm : mdeg

  w0 = 1 + 2/(13*(mdeg^2))
  temp1 = w0^2 - 1
  temp2 = sqrt(temp1)
  arg   = mdeg*log(w0 + temp2)
  w1    = (sinh(arg)*temp1) / (cosh(arg)*mdeg*temp2 - w0*sinh(arg))
  b1    = 1/((2*w0)^2)
  b2    = b1

  # stage-1
  gprev2 = copy(uprev)
  μs     = w1*b1
  gprev  = uprev + dt*μs*fsalfirst
  th2  = zero(eltype(u))
  th1  = μs
  z1   = w0
  z2   = one(eltype(u))
  dz1  = one(eltype(u))
  dz2  = zero(eltype(u))
  d2z1 = zero(eltype(u))
  d2z2 = zero(eltype(u))

  # stage 2 - mdeg
  for iter in 2:mdeg
    z   = 2*w0*z1 - z2
    dz  = 2*w0*dz1 - dz2 + 2*z1
    d2z = 2*w0*d2z1 - d2z2 + 4*dz1
    b   = d2z/(dz^2)
    νs  = 1 - z1*b1
    μ   = (2*w0*b)/b1
    ν   = - b/b2
    μs  = μ*w1/w0
    #using u as temporary storage
    u   = f(gprev, p, t + dt*th1)
    integrator.destats.nf += 1
    u   = μ*gprev + ν*gprev2  + (1 - μ - ν)*uprev + dt*μs*(u - νs*fsalfirst)
    th  = μ*th1 + ν*th2 + μs*(1 - νs)
    if (iter < mdeg)
      gprev2 = gprev
      gprev  = u
      th2  = th1
      th1  = th
      b2   = b1
      b1   = b
      z2   = z1
      z1   = z
      dz2  = dz1
      dz1  = dz
      d2z2 = d2z1
      d2z1 = d2z
    end
  end
  # error estimate
  if integrator.opts.adaptive
    tmp = 0.8*(uprev - u) + 0.4*dt*(fsalfirst + gprev)
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.u = u
end

function initialize!(integrator, cache::RKCCache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::RKCCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack k, tmp, gprev2, gprev, atmp = cache
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
  # The the number of degree for Chebyshev polynomial
  #maxm = max(2,Int(floor(sqrt(integrator.opts.internalnorm(integrator.opts.reltol,t)/10eps(t)))))
  maxm = 50
  mdeg = 1 + Int(floor(sqrt(1.54*abs(dt)*integrator.eigen_est + 1)))
  mdeg = (mdeg > maxm) ? maxm : mdeg

  w0 = 1 + 2/(13*(mdeg^2))
  temp1 = w0^2 - 1
  temp2 = sqrt(temp1)
  arg   = mdeg*log(w0 + temp2)
  w1    = (sinh(arg)*temp1) / (cosh(arg)*mdeg*temp2 - w0*sinh(arg))
  b1    = 1/((2*w0)^2)
  b2    = b1

  # stage-1
  @.. gprev2 = uprev
  μs     = w1*b1
  @.. gprev  = uprev + dt*μs*fsalfirst
  th2  = zero(eltype(u))
  th1  = μs
  z1   = w0
  z2   = one(eltype(u))
  dz1  = one(eltype(u))
  dz2  = zero(eltype(u))
  d2z1 = zero(eltype(u))
  d2z2 = zero(eltype(u))

  # stage 2 - mdeg
  for iter in 2:mdeg
    z   = 2*w0*z1 - z2
    dz  = 2*w0*dz1 - dz2 + 2*z1
    d2z = 2*w0*d2z1 - d2z2 + 4*dz1
    b   = d2z/(dz^2)
    νs  = 1 - z1*b1
    μ   = (2*w0*b)/b1
    ν   = - b/b2
    μs  = μ*w1/w0
    f(k, gprev, p, t + dt*th1)
    integrator.destats.nf += 1
    @.. u   = μ*gprev + ν*gprev2  + (1 - μ - ν)*uprev + dt*μs*(k - νs*fsalfirst)
    th  = μ*th1 + ν*th2 + μs*(1 - νs)
    if (iter < mdeg)
      gprev2 = gprev
      gprev  = u
      th2  = th1
      th1  = th
      b2   = b1
      b1   = b
      z2   = z1
      z1   = z
      dz2  = dz1
      dz1  = dz
      d2z2 = d2z1
      d2z1 = d2z
    end
  end
  # error estimate
  if integrator.opts.adaptive
    @.. tmp = 0.8*(uprev - u) + 0.4*dt*(fsalfirst + gprev)
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  f(integrator.fsallast, u, p, t+dt)
  integrator.destats.nf += 1
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::IRKCConstantCache)
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

function perform_step!(integrator,cache::IRKCConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p,alg,fsalfirst = integrator
  @unpack minm,du₁,du₂,nlsolver = cache
  @unpack f1, f2 = integrator.f
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

  # The the number of degree for Chebyshev polynomial
  #maxm = max(2,Int(floor(sqrt(integrator.opts.internalnorm(integrator.opts.reltol,t)/(10 *eps(integrator.opts.internalnorm(uprev,t)))))))
  maxm = 50
  mdeg = 1 + floor(Int, sqrt(1.54*abs(dt)*integrator.eigen_est + 1))
  mdeg = min(maxm, max(minm, mdeg))

  ω₀    = 1 + 2/(13 * (mdeg^2))
  temp₁ = ω₀^2 - 1
  temp₂ = sqrt(temp₁)
  θ     = mdeg*log(ω₀ + temp₂)
  ω₁    = (sinh(θ)*temp₁)/(cosh(θ)*mdeg*temp₂ - ω₀*sinh(θ))
  Bⱼ₋₂  = 1/(4 * ω₀^2)
  Bⱼ₋₁  = 1/ω₀

  #stage-1
  f1ⱼ₋₂  = du₁
  gprev2 = copy(uprev)
  μs     = ω₁*Bⱼ₋₁
  μs₁    = μs

  # initial guess for implicit part
  # if alg.extrapolant == :linear
  #   nlsolver.z = dt*du₁
  # else # :constant
  #   nlsolver.z = zero(u)
  # end

  nlsolver.z = dt*du₁

  nlsolver.tmp = uprev + dt*μs₁*du₂
  nlsolver.γ   = μs₁
  nlsolver.c   = μs
  markfirststage!(nlsolver)
  z = nlsolve!(nlsolver, integrator, cache, false)
  # nlsolvefail(nlsolver) && return
  gprev = nlsolver.tmp + μs₁*z

  Cⱼ₋₂   = zero(eltype(u))
  Cⱼ₋₁   = μs
  Tⱼ₋₁   = ω₀
  Tⱼ₋₂   = one(eltype(u))
  Tⱼ₋₁′  = one(eltype(u))
  Tⱼ₋₂′  = zero(eltype(u))
  Tⱼ₋₁″  = zero(eltype(u))
  Tⱼ₋₂″  = zero(eltype(u))

  #stage- 2...mdeg
  for iter in 2:mdeg
    Tⱼ   = 2*ω₀*Tⱼ₋₁ - Tⱼ₋₂
    Tⱼ′  = 2*ω₀*Tⱼ₋₁′ + 2*Tⱼ₋₁ - Tⱼ₋₂′
    Tⱼ″  = 2*ω₀*Tⱼ₋₁″ + 4*Tⱼ₋₁′ - Tⱼ₋₂″
    Bⱼ   = Tⱼ″/(Tⱼ′^2)
    μ    = (2*ω₀*Bⱼ)/Bⱼ₋₁
    ν    = - Bⱼ/Bⱼ₋₂
    μs   = (μ*ω₁)/ω₀
    νs   = -(1 - Tⱼ₋₁*Bⱼ₋₁)*μs
    Cⱼ   = μ*Cⱼ₋₁ + ν*Cⱼ₋₂ + μs + νs

    f1ⱼ₋₁  = f1(gprev, p, t+Cⱼ₋₁*dt)
    f2ⱼ₋₁  = f2(gprev, p, t+Cⱼ₋₁*dt)
    integrator.destats.nf += 1
    integrator.destats.nf2 += 1
    nlsolver.tmp = (1-μ-ν)*uprev + μ*gprev + ν*gprev2 + dt*μs*f2ⱼ₋₁ + dt*νs*du₂ + (νs - (1 -μ-ν)*μs₁)*dt*du₁ - ν*μs₁*dt*f1ⱼ₋₂
    nlsolver.z   = dt*f1ⱼ₋₁
    nlsolver.c   = Cⱼ
    z = nlsolve!(nlsolver, integrator, cache, false)
    # ignoring newton method's convergence failure
    # nlsolvefail(nlsolver) && return
    u = nlsolver.tmp + μs₁*z
    if (iter < mdeg)
      f1ⱼ₋₂= f1ⱼ₋₁
      gprev2 = gprev
      gprev  = u
      Cⱼ₋₂   = Cⱼ₋₁
      Cⱼ₋₁   = Cⱼ
      Bⱼ₋₂   = Bⱼ₋₁
      Bⱼ₋₁   = Bⱼ
      Tⱼ₋₂   = Tⱼ₋₁
      Tⱼ₋₁   = Tⱼ
      Tⱼ₋₂′  = Tⱼ₋₁′
      Tⱼ₋₁′  = Tⱼ′
      Tⱼ₋₂″  = Tⱼ₋₁″
      Tⱼ₋₁″  = Tⱼ″
    end
  end

  cache.du₁ = f1(u, p, t+dt)
  cache.du₂ = f2(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.destats.nf2 += 1
  # error estimate
  if isnewton(nlsolver) && integrator.opts.adaptive
    update_W!(integrator, cache, dt, false)
    tmp = dt*(0.5*(cache.du₂ - du₂) + (0.5 - μs₁)*(cache.du₁ - du₁))
    tmp = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  integrator.fsallast = cache.du₁ + cache.du₂
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::IRKCCache)
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

function perform_step!(integrator, cache::IRKCCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,alg = integrator
  @unpack gprev,gprev2,f1ⱼ₋₁,f1ⱼ₋₂,f2ⱼ₋₁,du₁,du₂,atmp,nlsolver = cache
  @unpack tmp,z = nlsolver
  @unpack minm = cache.constantcache
  @unpack f1, f2 = integrator.f

  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
  # The the number of degree for Chebyshev polynomial
  #maxm = max(2,int(floor(sqrt(integrator.opts.internalnorm(integrator.opts.reltol,t)/(10 *eps(integrator.opts.internalnorm(uprev,t)))))))
  maxm = 50
  mdeg = 1 + Int(floor(sqrt(1.54*abs(dt)*integrator.eigen_est + 1)))
  mdeg = (mdeg < minm) ? minm : mdeg
  mdeg = (mdeg >= maxm) ? maxm : mdeg

  ω₀    = 1 + 2/(13 *(mdeg^2))
  temp₁ = ω₀^2 - 1
  temp₂ = sqrt(temp₁)
  θ     = mdeg*log(ω₀ + temp₂)
  ω₁    = (sinh(θ)*temp₁)/(cosh(θ)*mdeg*temp₂ - ω₀*sinh(θ))
  Bⱼ₋₂  = 1 / (4 * ω₀^2)
  Bⱼ₋₁  = 1 / ω₀

  #stage-1
  f1ⱼ₋₂  = du₁
  @.. gprev2 = uprev
  μs     = ω₁*Bⱼ₋₁
  μs₁    = μs

  # initial guess
  # if alg.extrapolant == :linear
  #   @.. z = dt*du₁
  # else # :constant
  #   @.. z = zero(eltype(u))
  # end
  @.. nlsolver.z = dt*du₁

  @.. nlsolver.tmp = uprev + dt*μs₁*du₂
  nlsolver.γ   = μs₁
  nlsolver.c   = μs
  markfirststage!(nlsolver)
  z = nlsolve!(nlsolver, integrator, cache, false)
  # ignoring newton method's convergence failure
  # nlsolvefail(nlsolver) && return
  @.. gprev = nlsolver.tmp + μs₁*nlsolver.z

  Cⱼ₋₂   = zero(eltype(u))
  Cⱼ₋₁   = μs
  Tⱼ₋₁   = ω₀
  Tⱼ₋₂   = one(eltype(u))
  Tⱼ₋₁′  = one(eltype(u))
  Tⱼ₋₂′  = zero(eltype(u))
  Tⱼ₋₁″  = zero(eltype(u))
  Tⱼ₋₂″  = zero(eltype(u))

  #stage- 2...mdeg
  for iter in 2:mdeg
    Tⱼ   = 2*ω₀*Tⱼ₋₁ - Tⱼ₋₂
    Tⱼ′  = 2*ω₀*Tⱼ₋₁′ + 2*Tⱼ₋₁ - Tⱼ₋₂′
    Tⱼ″  = 2*ω₀*Tⱼ₋₁″ + 4*Tⱼ₋₁′ - Tⱼ₋₂″
    Bⱼ   = Tⱼ″/(Tⱼ′^2)
    μ    = (2*ω₀*Bⱼ)/Bⱼ₋₁
    ν    = - Bⱼ/Bⱼ₋₂
    μs   = (μ*ω₁)/ω₀
    νs   = -(1 - Tⱼ₋₁*Bⱼ₋₁)*μs
    Cⱼ   = μ*Cⱼ₋₁ + ν*Cⱼ₋₂ + μs + νs

    f1(f1ⱼ₋₁, gprev, p, t+Cⱼ₋₁*dt)
    f2(f2ⱼ₋₁, gprev, p, t+Cⱼ₋₁*dt)
    integrator.destats.nf += 1
    integrator.destats.nf2 += 1
    @.. nlsolver.tmp = (1-μ-ν)*uprev + μ*gprev + ν*gprev2 + dt*μs*f2ⱼ₋₁ + dt*νs*du₂ + (νs - (1-μ-ν)*μs₁)*dt*du₁ - ν*μs₁*dt*f1ⱼ₋₂
    @.. nlsolver.z   = dt*f1ⱼ₋₁
    nlsolver.c = Cⱼ

    z = nlsolve!(nlsolver, integrator, cache, false)
    # nlsolvefail(nlsolver) && return
    @.. u = nlsolver.tmp + μs₁*nlsolver.z
    if (iter < mdeg)
      @.. f1ⱼ₋₂  = f1ⱼ₋₁
      @.. gprev2 = gprev
      @.. gprev  = u
      Cⱼ₋₂   = Cⱼ₋₁
      Cⱼ₋₁   = Cⱼ
      Bⱼ₋₂   = Bⱼ₋₁
      Bⱼ₋₁   = Bⱼ
      Tⱼ₋₂   = Tⱼ₋₁
      Tⱼ₋₁   = Tⱼ
      Tⱼ₋₂′  = Tⱼ₋₁′
      Tⱼ₋₁′  = Tⱼ′
      Tⱼ₋₂″  = Tⱼ₋₁″
      Tⱼ₋₁″  = Tⱼ″
    end
  end

  @.. f1ⱼ₋₁ = du₁
  @.. f2ⱼ₋₁ = du₂
  f1(du₁, u, p, t+dt)
  f2(du₂, u, p, t+dt)
  integrator.destats.nf += 1
  integrator.destats.nf2 += 1
  # error estimate
  if isnewton(nlsolver) && integrator.opts.adaptive
    update_W!(integrator, cache, dt, false)
    @.. gprev = dt*0.5*(du₂ - f2ⱼ₋₁) + dt*(0.5 - μs₁)*(du₁ - f1ⱼ₋₁)
    nlsolver.cache.linsolve(vec(tmp),get_W(nlsolver),vec(gprev),false)
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  @.. integrator.fsallast = du₁ + du₂
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::ESERK4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ESERK4ConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack ms, Cᵤ, Cₑ= cache
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

  mdeg = Int(floor(sqrt(abs(dt)*integrator.eigen_est))+1)
  mdeg = (mdeg > 4000) ? 4000 : mdeg
  cache.mdeg = mdeg
  choosedeg_SERK!(integrator,cache)
  mdeg = cache.mdeg
  start = cache.start
  internal_deg = cache.internal_deg
  α = 2.0/(mdeg^2)

  u = zero(uprev)
  tmp = zero(uprev)

  for i in 1:4
    hᵢ = dt/i
    tᵢ = t
    Sᵢ = zero(u)
    uᵢ₋₁ = uprev
    uᵢ₋₂ = zero(u)
    for j in 1:i
      r  = tᵢ
      Sᵢ = (cache.Bᵢ[start])*uᵢ₋₁
      for st in 1:mdeg
        k = f(uᵢ₋₁, p, r)
        integrator.destats.nf += 1

        if st%internal_deg == 1
          uᵢ = uᵢ₋₁ + α*hᵢ*k
        else
          uᵢ = 2*uᵢ₋₁ - uᵢ₋₂ + 2*α*hᵢ*k
        end
        q = convert(Int, floor(st/internal_deg))
        r = tᵢ + α*(st^2 + q*internal_deg^2)*hᵢ
        Sᵢ = Sᵢ + (cache.Bᵢ[start+st])*uᵢ
        if st < mdeg
          uᵢ₋₂ = uᵢ₋₁
          uᵢ₋₁ = uᵢ
        end
      end

      if j < i
        tᵢ = tᵢ + hᵢ
        uᵢ₋₁ = Sᵢ
      end
    end

    u = u + Cᵤ[i]*Sᵢ
    integrator.opts.adaptive && (tmp = tmp + Cₑ[i]*Sᵢ)
  end

  u = u/6
  if integrator.opts.adaptive
    tmp = tmp/6
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.u = u
end

function initialize!(integrator, cache::ESERK4Cache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::ESERK4Cache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack uᵢ, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, k = cache
  @unpack ms, Cᵤ, Cₑ = cache.constantcache
  ccache = cache.constantcache
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

  mdeg = Int(floor(sqrt(abs(dt)*integrator.eigen_est))+1)
  mdeg = (mdeg > 4000) ? 4000 : mdeg
  ccache.mdeg = mdeg
  choosedeg_SERK!(integrator,cache)
  mdeg = ccache.mdeg
  start = ccache.start
  internal_deg = ccache.internal_deg
  α = 2.0/(mdeg^2)

  @.. u = zero(uprev)
  @.. tmp = zero(uprev)
  for i in 1:4
    hᵢ = dt/i
    tᵢ = t
    @.. Sᵢ = zero(u)
    @.. uᵢ₋₁ = uprev
    @.. uᵢ₋₂ = zero(u)
    for j in 1:i
      r  = tᵢ
      @.. Sᵢ = (cache.constantcache.Bᵢ[start])*uᵢ₋₁
      for st in 1:mdeg
        f(k, uᵢ₋₁, p, r)
        integrator.destats.nf += 1

        if st%internal_deg == 1
          @.. uᵢ = uᵢ₋₁ + α*hᵢ*k
        else
          @.. uᵢ = 2*uᵢ₋₁ - uᵢ₋₂ + 2*α*hᵢ*k
        end
        q = convert(Int, floor(st/internal_deg))
        r = tᵢ + α*(st^2 + q*internal_deg^2)*hᵢ
        @.. Sᵢ = Sᵢ + (cache.constantcache.Bᵢ[start+st])*uᵢ
        if st < mdeg
          @.. uᵢ₋₂ = uᵢ₋₁
          @.. uᵢ₋₁ = uᵢ
        end
      end

      if j < i
        tᵢ = tᵢ + hᵢ
        @.. uᵢ₋₁ = Sᵢ
      end
    end

    @.. u = u + Cᵤ[i]*Sᵢ
    integrator.opts.adaptive && (@.. tmp = tmp + Cₑ[i]*Sᵢ)
  end

  @.. u = u/6


  if integrator.opts.adaptive
    @.. tmp = tmp/6
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  f(integrator.fsallast, u, p, t+dt)
  integrator.destats.nf += 1
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::ESERK5ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ESERK5ConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack ms, Cᵤ, Cₑ, Bᵢ= cache
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

  mdeg = Int(floor(sqrt(abs(dt)*integrator.eigen_est/0.98))+1)
  mdeg = (mdeg > 2000) ? 2000 : mdeg
  cache.mdeg = mdeg
  choosedeg_SERK!(integrator,cache)
  mdeg = cache.mdeg
  start = cache.start
  internal_deg = cache.internal_deg
  α = 100.0/(49.0*mdeg^2)

  u = zero(uprev)
  tmp = zero(uprev)
  for i in 1:5
    hᵢ = dt/i
    tᵢ = t
    Sᵢ = zero(u)
    uᵢ₋₁ = uprev
    uᵢ₋₂ = zero(u)
    for j in 1:i
      r  = tᵢ
      Sᵢ = (Bᵢ[start])*uᵢ₋₁
      for st in 1:mdeg
        k = f(uᵢ₋₁, p, r)
        integrator.destats.nf += 1

        if st%internal_deg == 1
          uᵢ = uᵢ₋₁ + α*hᵢ*k
        else
          uᵢ = 2*uᵢ₋₁ - uᵢ₋₂ + 2*α*hᵢ*k
        end
        q = convert(Int, floor(st/internal_deg))
        r = tᵢ + α*(st^2 + q*internal_deg^2)*hᵢ
        Sᵢ = Sᵢ + (Bᵢ[start+st])*uᵢ
        if st < mdeg
          uᵢ₋₂ = uᵢ₋₁
          uᵢ₋₁ = uᵢ
        end
      end

      if j < i
        tᵢ = tᵢ + hᵢ
        uᵢ₋₁ = Sᵢ
      end
    end

    u = u + Cᵤ[i]*Sᵢ
    integrator.opts.adaptive && (tmp = tmp + Cₑ[i]*Sᵢ)
  end

  u = u/24
  if integrator.opts.adaptive
    tmp = tmp/24
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.u = u
end

function initialize!(integrator, cache::ESERK5Cache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::ESERK5Cache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack uᵢ, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, k = cache
  @unpack ms, Cᵤ, Cₑ, Bᵢ = cache.constantcache
  ccache = cache.constantcache
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

  mdeg = Int(floor(sqrt(abs(dt)*integrator.eigen_est/0.98))+1)
  mdeg = (mdeg > 2000) ? 2000 : mdeg
  ccache.mdeg = mdeg
  choosedeg_SERK!(integrator,cache)
  mdeg = ccache.mdeg
  start = ccache.start
  internal_deg = ccache.internal_deg
  α = 100.0/(49.0*mdeg^2)

  @.. u = zero(uprev)
  @.. tmp = zero(uprev)
  for i in 1:5
    hᵢ = dt/i
    tᵢ = t
    @.. Sᵢ = zero(u)
    @.. uᵢ₋₁ = uprev
    @.. uᵢ₋₂ = zero(u)
    for j in 1:i
      r  = tᵢ
      @.. Sᵢ = (Bᵢ[start])*uᵢ₋₁
      for st in 1:mdeg
        f(k, uᵢ₋₁, p, r)
        integrator.destats.nf += 1

        if st%internal_deg == 1
          @.. uᵢ = uᵢ₋₁ + α*hᵢ*k
        else
          @.. uᵢ = 2*uᵢ₋₁ - uᵢ₋₂ + 2*α*hᵢ*k
        end
        q = convert(Int, floor(st/internal_deg))
        r = tᵢ + α*(st^2 + q*internal_deg^2)*hᵢ
        @.. Sᵢ = Sᵢ + (Bᵢ[start+st])*uᵢ
        if st < mdeg
          @.. uᵢ₋₂ = uᵢ₋₁
          @.. uᵢ₋₁ = uᵢ
        end
      end

      if j < i
        tᵢ = tᵢ + hᵢ
        @.. uᵢ₋₁ = Sᵢ
      end
    end

    @.. u = u + Cᵤ[i]*Sᵢ
    integrator.opts.adaptive && (@.. tmp = tmp + Cₑ[i]*Sᵢ)
  end

  @.. u = u/24


  if integrator.opts.adaptive
    @.. tmp = tmp/24
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  f(integrator.fsallast, u, p, t+dt)
  integrator.destats.nf += 1
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::SERK2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::SERK2ConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack ms, Bᵢ= cache
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

  mdeg = Int(floor(sqrt(abs(dt)*integrator.eigen_est/0.8))+1)
  mdeg = (mdeg > 250) ? 250 : mdeg
  cache.mdeg = mdeg
  choosedeg_SERK!(integrator,cache)
  mdeg = cache.mdeg
  start = cache.start
  internal_deg = cache.internal_deg
  α = 1.0/(0.4*mdeg^2)

  uᵢ₋₁ = uprev
  uᵢ₋₂ = uprev
  Sᵢ   = Bᵢ[start]*uprev
  for i in 1:10
    k = f(uᵢ₋₁, p, t+(1+(i-1)*internal_deg^2)*α*dt)
    integrator.destats.nf += 1
    u    = uᵢ₋₁ + α*dt*k
    Sᵢ   = Sᵢ + Bᵢ[start + (i-1)*internal_deg + 1]*u
    uᵢ₋₂ = uᵢ₋₁
    uᵢ₋₁ = u
    for j in 2:internal_deg
      k = f(uᵢ₋₁, p, t+(j^2+(i-1)*internal_deg^2)*α*dt)
      integrator.destats.nf += 1
      u = 2*uᵢ₋₁ - uᵢ₋₂ + 2*α*dt*k
      Sᵢ= Sᵢ + Bᵢ[start+j+(i-1)*internal_deg]*u
      if j*i < mdeg
        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = u
      end
    end
  end
  u = Sᵢ
  k = f(u,p,t+dt)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    tmp = u - uprev - dt*k
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast = k
  integrator.u = u
end

function initialize!(integrator, cache::SERK2Cache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::SERK2Cache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, k = cache
  @unpack ms, Bᵢ = cache.constantcache
  ccache = cache.constantcache
  alg = unwrap_alg(integrator, true)
  alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

  mdeg = Int(floor(sqrt(abs(dt)*integrator.eigen_est/0.8))+1)
  mdeg = (mdeg > 250) ? 250 : mdeg
  ccache.mdeg = mdeg
  choosedeg_SERK!(integrator,cache)
  mdeg = ccache.mdeg
  start = ccache.start
  internal_deg = ccache.internal_deg
  α = 1.0/(0.4*mdeg^2)

  @.. uᵢ₋₁ = uprev
  @.. uᵢ₋₂ = uprev
  @.. Sᵢ   = Bᵢ[start]*uprev
  for i in 1:10
    f(k, uᵢ₋₁, p, t+(1+(i-1)*internal_deg^2)*α*dt)
    integrator.destats.nf += 1
    @.. u    = uᵢ₋₁ + α*dt*k
    @.. Sᵢ   = Sᵢ + Bᵢ[start + (i-1)*internal_deg + 1]*u
    @.. uᵢ₋₂ = uᵢ₋₁
    @.. uᵢ₋₁ = u
    for j in 2:internal_deg
      f(k, uᵢ₋₂, p, t+(j^2+(i-1)*internal_deg^2)*α*dt)
      integrator.destats.nf += 1
      @.. u = 2*uᵢ₋₁ - uᵢ₋₂ + 2*α*dt*k
      @.. Sᵢ= Sᵢ + Bᵢ[start+j+(i-1)*internal_deg]*u
      if j < mdeg
        @.. uᵢ₋₂ = uᵢ₋₁
        @.. uᵢ₋₁ = u
      end
    end
  end
  @.. u = Sᵢ
  f(k, u, p, t+dt)
  integrator.destats.nf += 1


  if integrator.opts.adaptive
    @.. tmp = u - uprev - dt*k
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast = k
  integrator.u = u
end
