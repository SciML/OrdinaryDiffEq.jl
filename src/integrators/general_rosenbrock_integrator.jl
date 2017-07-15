@inline function initialize!(integrator,cache::GenRosen4ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline function perform_step!(integrator,cache::GenRosen4ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack A,c,α,αEEst,stages = cache
  @unpack kk = cache
  # Calc First
  kk[1] = integrator.fsalfirst
  # Calc Middle
  for i = 2:stages-1
    utilde = zero(kk[1])
    for j = 1:i-1
      utilde = @muladd utilde + A[j,i]*kk[j]
    end
    kk[i] = f(@muladd(t+c[i]*dt),@muladd(uprev+dt*utilde));
  end
  #Calc Last
  utilde = zero(kk[1])
  for j = 1:stages-1
    utilde = @muladd utilde + A[j,end]*kk[j]
  end
  kk[end] = f(@muladd(t+c[end]*dt),@muladd(uprev+dt*utilde)); integrator.fsallast = kk[end] # Uses fsallast as temp even if not fsal
  # Accumulate Result
  utilde = α[1]*kk[1]
  for i = 2:stages
    utilde = @muladd utilde + α[i]*kk[i]
  end
  u = @muladd uprev + dt*utilde
  if integrator.opts.adaptive
    uEEst = αEEst[1]*kk[1]
    for i = 2:stages
      uEEst = @muladd uEEst + αEEst[i]*kk[i]
    end
    integrator.EEst = integrator.opts.internalnorm( dt*(utilde-uEEst)/@muladd(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
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
