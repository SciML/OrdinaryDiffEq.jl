function initialize!(integrator, cache::GenRosen4ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::GenRosen4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack A,c,α,αEEst,stages = cache
  @unpack kk = cache

  # Calc First
  kk[1] = integrator.fsalfirst

  # Calc Middle
  for i = 2:stages-1
    utilde = zero(kk[1])
    for j = 1:i-1
      utilde = utilde + A[j,i]*kk[j]
    end
    kk[i] = f(t+c[i]*dt, @.. uprev+dt*utilde)
    integrator.destats.nf += 1
  end

  # Calc Last
  utilde = zero(kk[1])
  for j = 1:stages-1
    utilde = utilde + A[j,end]*kk[j]
  end
  kk[end] = f(t+c[end]*dt, @.. uprev+dt*utilde)
  integrator.destats.nf += 1
  integrator.fsallast = kk[end] # Uses fsallast as temp even if not fsal

  # Accumulate Result
  utilde = α[1]*kk[1]
  for i = 2:stages
    utilde = utilde + α[i]*kk[i]
  end
  u = uprev + dt*utilde

  if integrator.opts.adaptive
    utilde = (α[1]-αEEst[1])*kk[1]
    for i = 2:stages
      utilde = utilde + (α[i]-αEEst[i])*kk[i]
    end
    utilde = dt*utilde
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  if !isfsal(integrator.alg.tableau)
    integrator.fsallast = f(u, p, t+dt)
    integrator.destats.nf += 1
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end
