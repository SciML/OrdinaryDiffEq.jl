@inline function initialize!(integrator,cache::ExplicitRKConstantCache)
  if isfsal(integrator.alg) # pre-start FSAL
    integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev)
  end
end

function perform_step!(integrator::ODEIntegrator,cache::ExplicitRKConstantCache)
  @unpack t,dt,uprev,u,f,k = integrator
  @unpack A,c,α,αEEst,stages = integrator.cache
  @unpack kk = integrator.cache
  # Calc First
  if isfsal(integrator.alg)
    kk[1] = integrator.fsalfirst
  else
    kk[1] = f(t,uprev)
  end
  # Calc Middle
  for i = 2:stages-1
    utilde = zero(kk[1])
    for j = 1:i-1
      utilde += A[j,i]*kk[j]
    end
    kk[i] = f(t+c[i]*dt,uprev+dt*utilde);
  end
  #Calc Last
  utilde = zero(kk[1])
  for j = 1:stages-1
    utilde += A[j,end]*kk[j]
  end
  kk[end] = f(t+c[end]*dt,uprev+dt*utilde); integrator.fsallast = kk[end] # Uses fsallast as temp even if not fsal
  # Accumulate Result
  utilde = α[1]*kk[1]
  for i = 2:stages
    utilde += α[i]*kk[i]
  end
  u = uprev + dt*utilde
  if integrator.opts.adaptive
    uEEst = αEEst[1]*kk[1]
    for i = 2:stages
      uEEst += αEEst[i]*kk[i]
    end
    integrator.EEst = abs( dt*(utilde-uEEst)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  if integrator.opts.calck
    k = kk[end]
  end
  @pack integrator = t,dt,u,k
end

function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O,algType<:ExplicitRK}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      perform_step!(integrator,integrator.cache)
      ode_loopfooter!(integrator)
      if isempty(integrator.tstops)
        break
      end
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

@inline function initialize!(integrator,cache::ExplicitRKCache)
  integrator.k = cache.kk[end]
  integrator.fsallast = cache.kk[end]
  integrator.fsalfirst = cache.kk[1]
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # Pre-start fsal
end

function perform_step!(integrator::ODEIntegrator,cache::ExplicitRKCache)
  @unpack t,dt,uprev,u,f,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack A,c,α,αEEst,stages = integrator.cache.tab
  @unpack kk,utilde,tmp,atmp,uEEst = integrator.cache
  # First
  if !isfsal(integrator.alg)
    f(t,uprev,kk[1])
  end
  # Middle
  for i = 2:stages-1
    for l in uidx
      utilde[l] = zero(kk[1][1])
    end
    for j = 1:i-1
      for l in uidx
        utilde[l] += A[j,i]*kk[j][l]
      end
    end
    for l in uidx
      tmp[l] = uprev[l]+dt*utilde[l]
    end
    f(t+c[i]*dt,tmp,kk[i])
  end
  #Last
  for l in uidx
    utilde[l] = zero(kk[1][1])
  end
  for j = 1:stages-1
    for l in uidx
      utilde[l] += A[j,end]*kk[j][l]
    end
  end
  for l in uidx
    u[l] = uprev[l]+dt*utilde[l]
  end
  f(t+c[end]*dt,u,kk[end]) #fsallast is tmp even if not fsal
  #Accumulate
  if !isfsal(integrator.alg)
    for i in uidx
      utilde[i] = α[1]*kk[1][i]
    end
    for i = 2:stages
      for l in uidx
        utilde[l] += α[i]*kk[i][l]
      end
    end
    for i in uidx
      u[i] = uprev[i] + dt*utilde[i]
    end
  end
  if integrator.opts.adaptive
    for i in uidx
      uEEst[i] = αEEst[1]*kk[1][i]
    end
    for i = 2:stages
      for j in uidx
        uEEst[j] += αEEst[i]*kk[i][j]
      end
    end
    for i in uidx
      atmp[i] = (dt*(utilde[i]-uEEst[i])/(integrator.opts.abstol+max(abs(uprev[i]),abs(u[i]))*integrator.opts.reltol))
    end
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  @pack integrator = t,dt,u,k
end

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O,algType<:ExplicitRK}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      perform_step!(integrator,integrator.cache)
      ode_loopfooter!(integrator)
      if isempty(integrator.tstops)
        break
      end
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end
