function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O,algType<:ExplicitRK}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  @unpack A,c,α,αEEst,stages = alg.tableau
  A = A' # Transpose A to column major looping
  kk = Array{ksEltype}(stages) # Not ks since that's for integrator.opts.dense
  local utilde::ksEltype
  if isfsal(integrator.alg) # pre-start FSAL
    fsalfirst = f(t,u)
  end
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      # Calc First
      if isfsal(integrator.alg)
        kk[1] = fsalfirst
      else
        kk[1] = f(t,u)
      end
      # Calc Middle
      for i = 2:stages-1
        utilde = zero(kk[1])
        for j = 1:i-1
          utilde += A[j,i]*kk[j]
        end
        kk[i] = f(t+c[i]*dt,u+dt*utilde);
      end
      #Calc Last
      utilde = zero(kk[1])
      for j = 1:stages-1
        utilde += A[j,end]*kk[j]
      end
      kk[end] = f(t+c[end]*dt,u+dt*utilde); fsallast = kk[end] # Uses fsallast as temp even if not fsal
      # Accumulate Result
      utilde = α[1]*kk[1]
      for i = 2:stages
        utilde += α[i]*kk[i]
      end
      utmp = u + dt*utilde
      if integrator.opts.adaptive
        uEEst = αEEst[1]*kk[1]
        for i = 2:stages
          uEEst += αEEst[i]*kk[i]
        end
        EEst = abs( dt*(utilde-uEEst)/(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol))
      end
      if integrator.opts.calck
        k = kk[end]
      end
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O,algType<:ExplicitRK}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  uidx = eachindex(u)
  @unpack A,c,α,αEEst,stages = alg.tableau
  A = A' # Transpose A to column major looping

  @unpack kk,utilde,tmp,atmp,uEEst = integrator.cache

  if integrator.opts.calck
    integrator.k = kk[end]
  end
  fsallast = kk[end]
  fsalfirst = kk[1]
  f(t,u,kk[1]) # pre-start fsal


  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      # First
      if !isfsal(integrator.alg)
        f(t,u,kk[1])
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
          tmp[l] = u[l]+dt*utilde[l]
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
        utmp[l] = u[l]+dt*utilde[l]
      end
      f(t+c[end]*dt,utmp,kk[end]) #fsallast is tmp even if not fsal
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
          utmp[i] = u[i] + dt*utilde[i]
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
          atmp[i] = (dt*(utilde[i]-uEEst[i])/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol))
        end
        EEst = integrator.opts.internalnorm(atmp)
      end
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end
