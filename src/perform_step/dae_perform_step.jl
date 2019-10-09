function initialize!(integrator, cache::DImplicitEulerCache) end
function initialize!(integrator, cache::DImplicitEulerConstantCache) end

@muladd function perform_step!(integrator, cache::DImplicitEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  alg = unwrap_alg(integrator, true)
  update_W!(integrator, cache, -inv(dt), repeat_step)
  @show cache.nlsolver.cache.W

  # initial guess
  guess = zero(u)

  # Wont work for non-mutable types
  zero_func! = (out, x) -> begin
                out = f((x .- uprev) ./ dt, x, p, t)
              end
  u = NLsolve.nlsolve(zero_func!, guess).zero

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
end


@muladd function perform_step!(integrator, cache::DImplicitEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack atmp,tmp = cache
  alg = unwrap_alg(integrator, true)
  update_W!(integrator, cache, -inv(dt), repeat_step)
  @show cache.nlsolver.cache.W

  # initial guess
  guess = tmp
  guess .= zero(eltype(u))

  zero_func! = (out, x) -> begin
                f(out, (x .- uprev) ./ dt, x, p, t)
              end
  u .= NLsolve.nlsolve(zero_func!, guess).zero

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
end
