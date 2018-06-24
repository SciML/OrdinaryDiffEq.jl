function initialize!(integrator, cache::ROCK2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ROCK2ConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  @unpack ms, fp1, fp2, recf = cache
  # WIP
  # The number of stage.
  mdeg = Int(floor(sqrt((1.5 + dt * integrator.eigen_est)/0.811) + 1))
  if mdeg >= 200
    dt = 0.8 * (200 ^ 2 * 0.811 - 1.5)/integrator.eigen_est
    mdeg = 200
  end
  cache.mdeg = max(mdeg, 3) - 2
  cache.mdeg != cache.mdegprev && choosedeg!(cache)
  err = 0
  # recurrence
  # for the first stage
  temp1 = dt * recf[cache.recind]
  ci1 = t + temp1
  ci2 = t + temp1
  ci3 = t
  gprev2 = uprev
  gprev = uprev + temp1 * fsalfirst
  ms[cache.mdeg] < 2 && ( u = gprev )
  # for the second to the ms[cache.mdeg] th stages
  for i in 2:ms[cache.mdeg]
    temp1 = dt * recf[cache.recind + 2 * (i - 2) + 1]
    temp3 = -recf[cache.recind + 2 * (i - 2) + 2]
    temp2 = 1 - temp3
    ci1 = temp1 + temp2 * ci2 + temp3 * ci3
    u = temp1 * u + temp2 * gprev + temp3 * gprev2
    i < ms[cache.mdeg] && (gprev2 = gprev; gprev = u)
    ci3 = ci2
    ci2 = ci1
  end # end if
  # two-stage finishing procedure.
  temp1 = dt * fp1[cache.mdeg]
  temp2 = dt * fp2[cache.mdeg]
  gprev = u + temp1 * gprev2
  ci1 += temp1
  # error estimate
  if integrator.opts.adaptive
    temp3 = temp2 * (u - gprev2)
    u = gprev + temp1 * u + temp3
    calculate_residuals(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  #
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end
