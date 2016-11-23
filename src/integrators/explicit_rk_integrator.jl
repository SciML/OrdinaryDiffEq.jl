function ode_solve{uType<:Number,uEltype,N,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5}(integrator::ODEIntegrator{ExplicitRK,uType,uEltype,N,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5})
  @ode_preamble
  local A::Matrix{uEltypeNoUnits}
  local c::Vector{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local αEEst::Vector{uEltypeNoUnits}
  local stages::Int
  @unpack A,c,α,αEEst,stages,fsal = integrator.tableau
  A = A' # Transpose A to column major looping
  kk = Array{rateType}(stages) # Not ks since that's for dense
  local utilde::rateType
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      # Calc First
      if fsal
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
      if adaptive
        utmp = u + dt*utilde
        uEEst = αEEst[1]*kk[1]
        for i = 2:stages
          uEEst += αEEst[i]*kk[i]
        end
        EEst = abs( dt*(utilde-uEEst)/(abstol+max(abs(u),abs(utmp))*reltol))
      else
        u = u + dt*utilde
      end
      if calck
        k = kk[end]
      end
      @ode_loopfooter
    end
  end
  @ode_postamble
end

function ode_solve{uType<:AbstractArray,uEltype,N,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5}(integrator::ODEIntegrator{ExplicitRK,uType,uEltype,N,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5})
  @ode_preamble
  local A::Matrix{uEltypeNoUnits}
  local c::Vector{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local αEEst::Vector{uEltypeNoUnits}
  local stages::Int
  uidx = eachindex(u)
  @unpack A,c,α,αEEst,stages,fsal = integrator.tableau
  A = A' # Transpose A to column major looping
  kk = Vector{rateType}(0)
  for i = 1:stages
    push!(kk,similar(rate_prototype))
  end
  utilde = similar(rate_prototype)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  utmp = zeros(u)
  uEEst = similar(rate_prototype)
  fsallast = kk[end]
  fsalfirst = kk[1]
  if calck
    k = kk[end]
  end
  f(t,u,kk[1]) # pre-start fsal
  if custom_callback
    cache = (u,tmp,utilde,uEEst,atmp,uprev,kprev,utmp,kk...)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      # First
      if !fsal
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
      if !fsal
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
      if adaptive
        for i in uidx
          uEEst[i] = αEEst[1]*kk[1][i]
        end
        for i = 2:stages
          for j in uidx
            uEEst[j] += αEEst[i]*kk[i][j]
          end
        end
        for i in uidx
          atmp[i] = (dt*(utilde[i]-uEEst[i])/(abstol+max(abs(u[i]),abs(utmp[i]))*reltol))
        end
        EEst = internalnorm(atmp)
      else
        recursivecopy!(u,utmp)
      end
      @ode_loopfooter
    end
  end
  @ode_postamble
end
