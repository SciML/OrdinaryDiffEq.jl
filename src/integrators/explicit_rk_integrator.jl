function initialize!(integrator, cache::ExplicitRKConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.t, integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ExplicitRKConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack A,c,α,αEEst,stages = cache
  @unpack kk = cache

  # Calc First
  kk[1] = integrator.fsalfirst

  # Calc Middle
  for i = 2:stages-1
    utilde = zero(kk[1])
    for j = 1:i-1
      utilde = @. utilde + A[j,i]*kk[j]
    end
    kk[i] = f(t+c[i]*dt, @. uprev + dt*utilde);
  end

  #Calc Last
  utilde = zero(kk[1])
  for j = 1:stages-1
    utilde = @. utilde + A[j,end]*kk[j]
  end
  kk[end] = f(t+c[end]*dt, @. uprev + dt*utilde)
  integrator.fsallast = kk[end] # Uses fsallast as temp even if not fsal

  # Accumulate Result
  utilde = α[1]*kk[1]
  for i = 2:stages
    utilde = @. utilde + α[i]*kk[i]
  end
  u = @. uprev + dt*utilde

  if integrator.opts.adaptive
    utilde = (α[1]-αEEst[1])*kk[1]
    for i = 2:stages
      utilde = @. utilde + (α[i]-αEEst[i])*kk[i]
    end
    utilde = dt*utilde
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  if !isfsal(integrator.alg.tableau)
    integrator.fsallast = f(t+dt, u)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::ExplicitRKCache)
  integrator.kshortsize = 2
  integrator.fsallast = cache.fsallast
  integrator.fsalfirst = cache.kk[1]
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.t, integrator.uprev, integrator.fsalfirst) # Pre-start fsal
end

@muladd function perform_step!(integrator, cache::ExplicitRKCache, repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack A,c,α,αEEst,stages = cache.tab
  @unpack kk,utilde,tmp,atmp = cache

  # Middle
  for i = 2:stages-1
    @. utilde = zero(kk[1][1])
    for j = 1:i-1
      @. utilde = utilde + A[j,i]*kk[j]
    end
    @. tmp = uprev+dt*utilde
    f(t+c[i]*dt, tmp, kk[i])
  end

  #Last
  @. utilde = zero(kk[1][1])
  for j = 1:stages-1
    @. utilde = utilde + A[j,end]*kk[j]
  end
  @. u = uprev + dt*utilde
  f(t+c[end]*dt, u, kk[end]) #fsallast is tmp even if not fsal

  #Accumulate
  if !isfsal(integrator.alg.tableau)
    @. utilde = α[1]*kk[1]
    for i = 2:stages
      @. utilde = utilde + α[i]*kk[i]
    end
    @. u = uprev + dt*utilde
  end

  if integrator.opts.adaptive
    @. utilde = (α[1]-αEEst[1])*kk[1]
    for i = 2:stages
      @. utilde = utilde + (α[i]-αEEst[i])*kk[i]
    end
    @. utilde = dt*utilde
    calculate_residuals!(atmp, utilde, uprev, u,
                         integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  if !isfsal(integrator.alg.tableau)
    f(t+dt, u, integrator.fsallast)
  end
end
