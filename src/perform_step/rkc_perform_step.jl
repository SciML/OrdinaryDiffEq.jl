function initialize!(integrator, cache::ROCK2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ROCK2ConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack ms, fp1, fp2, recf = cache
  maxeig!(integrator, cache)
  # The the number of degree for Chebyshev polynomial
  mdeg = Int(floor(sqrt((1.5 + dt * integrator.eigen_est)/0.811) + 1))
  mdeg = min(mdeg, 200)
  cache.mdeg = max(mdeg, 3) - 2
  cache.mdeg != cache.mdegprev && choosedeg!(cache)
  # recurrence
  # for the first stage
  temp1 = dt * recf[cache.recind]
  ci1 = t + temp1
  ci2 = t + temp1
  ci3 = t
  gprev2 = copy(uprev)
  gprev = uprev + temp1 * fsalfirst
  ms[cache.mdeg] < 2 && ( u = gprev )
  # for the second to the ms[cache.mdeg] th stages
  for i in 2:ms[cache.mdeg]
    μ, κ = recf[cache.recind + (i - 2)*2 + 1], recf[cache.recind + (i - 2)*2 + 2]
    ν = -1 - κ
    dtμ = dt*μ
    ci1 = dtμ - ν * ci2 - κ * ci3
    u = dtμ * u - ν * gprev - κ * gprev2
    i < ms[cache.mdeg] && (gprev2 = gprev; gprev = u)
    ci3 = ci2
    ci2 = ci1
  end # end if
  # two-stage finishing procedure.
  temp1 = dt * fp1[cache.mdeg]
  temp2 = dt * fp2[cache.mdeg]
  gprev2 = f(u, p, ci1)
  integrator.destats.nf += 1
  gprev = u + temp1 * gprev2
  ci1 += temp1
  u = f(gprev, p, ci1)
  integrator.destats.nf += 1
  temp3 = temp2 * (u - gprev2)
  u = gprev + temp1 * u + temp3
  # error estimate
  if integrator.opts.adaptive
    atmp = calculate_residuals(temp3, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
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
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::ROCK2Cache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack k, tmp, gprev2, gprev, atmp = cache
  @unpack ms, fp1, fp2, recf = cache.constantcache
  ccache = cache.constantcache
  maxeig!(integrator, cache)
  # The the number of degree for Chebyshev polynomial
  mdeg = Int(floor(sqrt((1.5 + dt * integrator.eigen_est)/0.811) + 1))
  mdeg = min(mdeg, 200)
  ccache.mdeg = max(mdeg, 3) - 2
  ccache.mdeg != ccache.mdegprev && choosedeg!(cache)
  # recurrence
  # for the first stage
  temp1 = dt * recf[ccache.recind]
  ci1 = t + temp1
  ci2 = t + temp1
  ci3 = t
  @.. gprev2 = uprev
  @.. gprev = uprev + temp1 * fsalfirst
  ms[ccache.mdeg] < 2 && ( @.. u = gprev )
  # for the second to the ms[ccache.mdeg] th stages
  for i in 2:ms[ccache.mdeg]
    μ, κ = recf[ccache.recind + (i - 2)*2 + 1], recf[ccache.recind + (i - 2)*2 + 2]
    ν = κ - 1
    temp1 = dt * μ
    temp2 = 1 + κ
    temp3 = -κ
    ci1 = temp1 + temp2 * ci2 + temp3 * ci3
    @.. u = temp1 * u + temp2 * gprev + temp3 * gprev2
    i < ms[ccache.mdeg] && (gprev2 .= gprev; gprev .= u)
    ci3 = ci2
    ci2 = ci1
  end # end if
  # two-stage finishing procedure.
  temp1 = dt * fp1[ccache.mdeg]
  temp2 = dt * fp2[ccache.mdeg]
  f(k, u, p, ci1)
  integrator.destats.nf += 1
  @.. gprev = u + temp1 * k
  @.. tmp = -temp2 * k
  ci1 += temp1
  f(k, gprev, p, ci1)
  integrator.destats.nf += 1
  @.. tmp += temp2 * k
  @.. u = gprev + temp1 * k + tmp
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
  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ROCK4ConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack ms, fpa, fpb, fpbe, recf = cache
  maxeig!(integrator, cache)
  # The the number of degree for Chebyshev polynomial
  mdeg = Int(floor(sqrt((3 + dt * integrator.eigen_est)/0.353) + 1))
  if mdeg >= 152
    mdeg = 152
  end
  cache.mdeg = max(mdeg, 5) - 4
  cache.mdeg != cache.mdegprev && choosedeg!(cache)
  # recurrence
  # for the first stage
  temp1 = dt * recf[cache.recind]
  ci1 = t + temp1
  ci2 = t + temp1
  ci3 = t
  gprev2 = copy(uprev)
  gprev = uprev + temp1 * fsalfirst
  ms[cache.mdeg] < 2 && ( u = gprev )
  # for the second to the ms[cache.mdeg] th stages
  for i in 2:ms[cache.mdeg]
    μ, κ = recf[cache.recind + (i - 2)*2 + 1], recf[cache.recind + (i - 2)*2 + 2]
    ν = -1 - κ
    dtμ = dt*μ
    ci1 = dtμ - ν * ci2 - κ * ci3
    u = dtμ * u - ν * gprev - κ * gprev2
    i < ms[cache.mdeg] && (gprev2 = gprev; gprev = u)
    ci3 = ci2
    ci2 = ci1
  end
  # 4-stage finishing procedure.
  # Stage-1
  temp1 = dt * fpa[cache.mdeg][1]
  gprev = f(u, p, ci1)
  integrator.destats.nf += 1
  gprev3 = u + temp1 * gprev
  # Stage-2
  ci2 = ci1 + temp1
  temp1 = dt * fpa[cache.mdeg][2];
  temp2 = dt * fpa[cache.mdeg][3];
  gprev2 = f(gprev3, p, ci1)
  integrator.destats.nf += 1
  gprev4 = u + temp1 * gprev + temp2 * gprev2
  # Stage-3
  ci2 = ci1 + temp1 +temp2
  temp1 = dt * fpa[cache.mdeg][4]
  temp2 = dt * fpa[cache.mdeg][5]
  temp3 = dt * fpa[cache.mdeg][6]
  gprev3 = f(gprev4, p, ci2)
  integrator.destats.nf += 1
  gprev5 = u + temp1 * gprev + temp2 * gprev2 + temp3 * gprev3
  #Stage-4
  ci2 = ci1 + temp1 + temp2 + temp3
  temp1 = dt * fpb[cache.mdeg][1]
  temp2 = dt * fpb[cache.mdeg][2]
  temp3 = dt * fpb[cache.mdeg][3]
  temp4 = dt * fpb[cache.mdeg][4]
  gprev4 = f(gprev5, p, ci2)
  integrator.destats.nf += 1
  u = u + temp1 * gprev + temp2 * gprev2 + temp3 * gprev3 + temp4 * gprev4
  #Error estimate (embedded method of order 3)
  temp1 = dt * fpbe[cache.mdeg][1] - temp1
  temp2 = dt * fpbe[cache.mdeg][2] - temp2
  temp3 = dt * fpbe[cache.mdeg][3] - temp3
  temp4 = dt * fpbe[cache.mdeg][4] - temp4
  temp5 = dt * fpbe[cache.mdeg][5]
  gprev5 = f(u, p, t + dt)
  integrator.destats.nf += 1
  temp5 = temp1 * gprev + temp2 * gprev2 + temp3 * gprev3 + temp4 * gprev4 + temp5 * gprev5
  if integrator.opts.adaptive
    atmp = calculate_residuals(temp5, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1
  integrator.u = u
end

function initialize!(integrator, cache::ROCK4Cache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::ROCK4Cache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack gprev, gprev2, gprev3, gprev4, gprev5, tmp, atmp, k, k2, k3, k4, k5 = cache
  @unpack ms, fpa, fpb, fpbe, recf = cache.constantcache
  ccache = cache.constantcache
  maxeig!(integrator, cache)
  # The the number of degree for Chebyshev polynomial
  mdeg = Int(floor(sqrt((3 + dt * integrator.eigen_est)/0.353) + 1))
  if mdeg >= 152
    mdeg = 152
  end
  ccache.mdeg = max(mdeg, 5) - 4
  ccache.mdeg != ccache.mdegprev && choosedeg!(cache)
  # recurrence
  # for the first stage
  temp1 = dt * recf[ccache.recind]
  ci1 = t + temp1
  ci2 = t + temp1
  ci3 = t
  @.. gprev2 = uprev
  @.. gprev = uprev + temp1 * fsalfirst
  ms[ccache.mdeg] < 2 && ( @.. u = gprev )
  # for the second to the ms[ccache.mdeg] th stages
  for i in 2:ms[ccache.mdeg]
    μ, κ = recf[ccache.recind + (i - 2)*2 + 1], recf[ccache.recind + (i - 2)*2 + 2]
    ν = κ - 1
    temp1 = dt * μ
    temp2 = 1 + κ
    temp3 = -κ
    ci1 = temp1 + temp2 * ci2 + temp3 * ci3
    @.. u = temp1 * u + temp2 * gprev + temp3 * gprev2
    i < ms[ccache.mdeg] && (gprev2 .= gprev; gprev .= u)
    ci3 = ci2
    ci2 = ci1
  end
  # 4-stage finishing procedure.
  # Stage-1
  temp1 = dt * fpa[ccache.mdeg][1]
  f(k, u, p, ci1)
  @.. gprev3 = u + temp1 * k
  # Stage-2
  ci2 = ci1 + temp1
  temp1 = dt * fpa[ccache.mdeg][2];
  temp2 = dt * fpa[ccache.mdeg][3];
  f(k2, gprev3, p, ci1)
  @.. gprev4 = u + temp1 * k + temp2 * k2
  # Stage-3
  ci2 = ci1 + temp1 +temp2
  temp1 = dt * fpa[ccache.mdeg][4]
  temp2 = dt * fpa[ccache.mdeg][5]
  temp3 = dt * fpa[ccache.mdeg][6]
  f(k3, gprev4, p, ci2)
  @.. gprev5 = u + temp1 * k + temp2 * k2 + temp3 * k3
  #Stage-4
  ci2 = ci1 + temp1 + temp2 + temp3
  temp1 = dt * fpb[ccache.mdeg][1]
  temp2 = dt * fpb[ccache.mdeg][2]
  temp3 = dt * fpb[ccache.mdeg][3]
  temp4 = dt * fpb[ccache.mdeg][4]
  f(k4, gprev5, p, ci2)
  @.. u = u + temp1 * k + temp2 * k2 + temp3 * k3 + temp4 * k4
  #Error estimate (embedded method of order 3)
  temp1 = dt * fpbe[ccache.mdeg][1] - temp1
  temp2 = dt * fpbe[ccache.mdeg][2] - temp2
  temp3 = dt * fpbe[ccache.mdeg][3] - temp3
  temp4 = dt * fpbe[ccache.mdeg][4] - temp4
  temp5 = dt * fpbe[ccache.mdeg][5]
  f(k5, u, p, t + dt)
  integrator.destats.nf += 5
  @.. tmp = temp1 * k + temp2 * k2 + temp3 * k3 + temp4 * k4 + temp5 * k5
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
  maxeig!(integrator, cache)
  # The the number of degree for Chebyshev polynomial
  maxm = max(2,Int(floor(sqrt(integrator.opts.internalnorm(integrator.opts.reltol,t)/(10*eps(integrator.opts.internalnorm(uprev,t)))))))
  mdeg = 1 + Int(floor(sqrt(1.54*dt*integrator.eigen_est + 1)))
  if mdeg >= maxm
    mdeg = maxm
  end

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
  maxeig!(integrator, cache)
  # The the number of degree for Chebyshev polynomial
  maxm = max(2,Int(floor(sqrt(integrator.opts.internalnorm(integrator.opts.reltol,t)/(10*eps(integrator.opts.internalnorm(uprev,t)))))))
  mdeg = 1 + Int(floor(sqrt(1.54*dt*integrator.eigen_est + 1)))
  if mdeg >= maxm
    mdeg = maxm
  end

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
  maxeig!(integrator, cache)

  # The the number of degree for Chebyshev polynomial
  maxm = max(2,Int(floor(sqrt(integrator.opts.internalnorm(integrator.opts.reltol,t)/(10 *eps(integrator.opts.internalnorm(uprev,t)))))))
  mdeg = 1 + Int(floor(sqrt(1.54*dt*integrator.eigen_est + 1)))
  mdeg = (mdeg < minm) ? minm : mdeg
  mdeg = (mdeg >= maxm) ? maxm : mdeg

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

  update_W!(integrator, cache, μs₁*dt, false)
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
  z = nlsolve!(integrator, cache)
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
    z = nlsolve!(integrator, cache)
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
  integrator.fsallast = cache.k
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
  @unpack tmp,gprev,gprev2,k,f1ⱼ₋₁,f1ⱼ₋₂,f2ⱼ₋₁,du₁,du₂,z,W,atmp,nlsolver = cache
  @unpack minm = cache.constantcache
  @unpack f1, f2 = integrator.f

  maxeig!(integrator, cache)
  # The the number of degree for Chebyshev polynomial
  maxm = max(2,Int(floor(sqrt(integrator.opts.internalnorm(integrator.opts.reltol,t)/(10 *eps(integrator.opts.internalnorm(uprev,t)))))))
  mdeg = 1 + Int(floor(sqrt(1.54*dt*integrator.eigen_est + 1)))
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

  update_W!(integrator, cache, μs₁*dt, false)
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
  z = nlsolve!(integrator, cache)
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

    z = nlsolve!(integrator, cache)
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
    if DiffEqBase.has_invW(f)
        mul!(vec(tmp),get_W(nlsolver),vec(gprev))
      else
        cache.linsolve(vec(tmp),get_W(nlsolver),vec(gprev),false)
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  @.. integrator.fsallast = du₁ + du₂
  integrator.k[1] = integrator.fsalfirst
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
  maxeig!(integrator, cache)

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
  maxeig!(integrator, cache)

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

function initialize!(integrator, cache::SERK2v2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::SERK2v2ConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack ms, Bᵢ= cache
  maxeig!(integrator, cache)

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

function initialize!(integrator, cache::SERK2v2Cache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::SERK2v2Cache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, k = cache
  @unpack ms, Bᵢ = cache.constantcache
  ccache = cache.constantcache
  maxeig!(integrator, cache)

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
