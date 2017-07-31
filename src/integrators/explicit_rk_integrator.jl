@inline function initialize!(integrator,cache::ExplicitRKConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline @muladd function perform_step!(integrator,cache::ExplicitRKConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
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
  kk[end] = f(t+c[end]*dt, @. uprev + dt*utilde); integrator.fsallast = kk[end] # Uses fsallast as temp even if not fsal
  # Accumulate Result
  utilde = α[1]*kk[1]
  for i = 2:stages
    utilde = @. utilde + α[i]*kk[i]
  end
  u = @. uprev + dt*utilde
  if integrator.opts.adaptive
    uEEst = αEEst[1]*kk[1]
    for i = 2:stages
      uEEst = @. uEEst + αEEst[i]*kk[i]
    end
    tmp = @. dt*(utilde-uEEst)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  if isfsal(integrator.alg.tableau)
    integrator.fsallast = kk[end]
  else
    integrator.fsallast = f(t+dt,u)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end

@inline function initialize!(integrator,cache::ExplicitRKCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.fsallast = cache.fsallast
  integrator.fsalfirst = cache.kk[1]
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # Pre-start fsal
end

#=
@inline @muladd function perform_step!(integrator,cache::ExplicitRKCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  @unpack A,c,α,αEEst,stages = cache.tab
  @unpack kk,utilde,tmp,atmp,uEEst = cache
  # Middle
  for i = 2:stages-1
    @. utilde = zero(kk[1][1])
    for j = 1:i-1
      @. utilde = utilde + A[j,i]*kk[j]
    end
    @. tmp = uprev+dt*utilde
    f(t+c[i]*dt,tmp,kk[i])
  end
  #Last
  @. utilde = zero(kk[1][1])
  for j = 1:stages-1
    @. utilde = utilde + A[j,end]*kk[j]
  end
  @. u = uprev+dt*utilde
  f(t+c[end]*dt),u,kk[end]) #fsallast is tmp even if not fsal
  #Accumulate
  if !isfsal(integrator.alg.tableau)
    @. utilde = α[1]*kk[1]
    for i = 2:stages
      @. utilde = utilde + α[i]*kk[i]
    end
    @. u = uprev + dt*utilde
  end
  if integrator.opts.adaptive
    @. uEEst = αEEst[1]*kk[1]
    for i = 2:stages
      @. uEEst = uEEst + αEEst[i]*kk[i]
    end
    @. atmp = (dt*(utilde-uEEst)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  if !isfsal(integrator.alg.tableau)
    f(t+dt,u,integrator.fsallast)
  end
  @pack integrator = t,dt,u
end
=#

@inline @muladd function perform_step!(integrator,cache::ExplicitRKCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  @unpack A,c,α,αEEst,stages = cache.tab
  @unpack kk,utilde,tmp,atmp,uEEst = cache
  # Middle
  for i = 2:stages-1
    @tight_loop_macros for l in uidx
      @inbounds utilde[l] = zero(kk[1][1])
    end
    for j = 1:i-1
      @tight_loop_macros for l in uidx
        @inbounds utilde[l] = utilde[l] + A[j,i]*kk[j][l]
      end
    end
    @tight_loop_macros for l in uidx
      @inbounds tmp[l] = uprev[l]+dt*utilde[l]
    end
    f(t+c[i]*dt,tmp,kk[i])
  end
  #Last
  @tight_loop_macros for l in uidx
    @inbounds utilde[l] = zero(kk[1][1])
  end
  for j = 1:stages-1
    @tight_loop_macros for l in uidx
      @inbounds utilde[l] = utilde[l] + A[j,end]*kk[j][l]
    end
  end
  @tight_loop_macros for l in uidx
    @inbounds u[l] = uprev[l]+dt*utilde[l]
  end
  f(t+c[end]*dt,u,kk[end]) #fsallast is tmp even if not fsal
  #Accumulate
  if !isfsal(integrator.alg.tableau)
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = α[1]*kk[1][i]
    end
    for i = 2:stages
      @tight_loop_macros for l in uidx
        @inbounds utilde[l] = utilde[l] + α[i]*kk[i][l]
      end
    end
    @tight_loop_macros for i in uidx
      @inbounds u[i] = uprev[i] + dt*utilde[i]
    end
  end
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds uEEst[i] = αEEst[1]*kk[1][i]
    end
    for i = 2:stages
      @tight_loop_macros for j in uidx
        @inbounds uEEst[j] = uEEst[j] + αEEst[i]*kk[i][j]
      end
    end
    @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds atmp[i] = dt*(utilde[i]-uEEst[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  if !isfsal(integrator.alg.tableau)
    f(t+dt,u,integrator.fsallast)
  end
  @pack integrator = t,dt,u
end
