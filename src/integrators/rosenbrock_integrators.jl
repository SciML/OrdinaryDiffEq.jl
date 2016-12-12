function ode_solve{uType<:AbstractArray,uEltype,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5}(integrator::ODEIntegrator{Rosenbrock23,uType,uEltype,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  local k₁::uType = similar(u)
  local k₂ = similar(u)
  local k₃::uType = similar(u)
  local tmp::uType
  const kshortsize = 1
  function vecf(t,u,du)
    f(t,reshape(u,sizeu...),reshape(du,sizeu...))
    u = vec(u)
    du = vec(du)
  end
  function vecfreturn(t,u,du)
    f(t,reshape(u,sizeu...),reshape(du,sizeu...))
    return vec(du)
  end
  du1 = zeros(u)
  du2 = zeros(u)
  f₀ = similar(u)
  f₁ = similar(u)
  vectmp3 = similar(vec(u))
  utmp = similar(u); vectmp2 = similar(vec(u))
  dT = similar(u); vectmp = similar(vec(u))
  J = zeros(uEltype,length(u),length(u))
  W = similar(J); tmp2 = similar(u)
  uidx = eachindex(u)
  jidx = eachindex(J)
  f(t,u,fsalfirst)
  if calck
    k = fsalfirst
  end
  cache = (u,du1,du2,uprev,kprev,f₀,f₁,vectmp3,utmp,vectmp2,dT,vectmp,tmp2,k)
  jaccache = (jidx,J,W)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      #if autodiff
        ForwardDiff.derivative!(dT,(t)->vecfreturn(t,u,du2),t) # Time derivative
        ForwardDiff.jacobian!(J,(du1,u)->vecf(t,u,du1),vec(du1),vec(u))
      #else
      #  Calculus.finite_difference!((t)->vecfreturn(t,u,du2),[t],dT)
      #  Calculus.finite_difference_jacobian!((du1,u)->vecf(t,u,du1),vec(u),vec(du1),J)
      #end

      W[:] = I-dt*d*J # Can an allocation be cut here?
      @into! vectmp = W\vec(fsalfirst + dt*d*dT)
      k₁ = reshape(vectmp,sizeu...)
      for i in uidx
        utmp[i]=u[i]+dt*k₁[i]/2
      end
      f(t+dt/2,utmp,f₁)
      @into! vectmp2 = W\vec(f₁-k₁)
      tmp = reshape(vectmp2,sizeu...)
      for i in uidx
        k₂[i] = tmp[i] + k₁[i]
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + dt*k₂[i]
        end
        f(t+dt,utmp,fsallast)
        @into! vectmp3 = W\vec(fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT)
        k₃ = reshape(vectmp3,sizeu...)
        for i in uidx
          tmp2[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)./(abstol+u[i]*reltol)
        end
        EEst = internalnorm(tmp2)
      else
        for i in uidx
          u[i] = u[i] + dt*k₂[i]
        end
        f(t,u,fsallast)
      end
      @ode_loopfooter
      recursivecopy!(fsalfirst,fsallast)
    end
  end
  @ode_postamble
end

function ode_solve{uType<:Number,uEltype,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5}(integrator::ODEIntegrator{Rosenbrock23,uType,uEltype,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  local dT::uType
  local J::uType
  local k₁::uType
  local f₁::uType
  local k₂::uType
  local k₃::uType
  const kshortsize = 1
  fsalfirst = f(t,u)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      # Time derivative
      dT = ForwardDiff.derivative((t)->f(t,u),t)
      J = ForwardDiff.derivative((u)->f(t,u),u)
      W = I-dt*d*J

      if calck
        k = fsalfirst
      end
      k₁ = W\(fsalfirst + dt*d*dT)
      f₁ = f(t+dt/2,u+dt*k₁/2)
      k₂ = W\(f₁-k₁) + k₁
      if adaptive
        utmp = u + dt*k₂
        fsallast = f(t+dt,utmp)
        k₃ = W\(fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT)
        EEst = abs((dt*(k₁ - 2k₂ + k₃)/6)./(abstol+u*reltol))
      else
        u = u + dt*k₂
        fsallast = f(t,u)
      end
      @ode_loopfooter
      fsalfirst = fsallast
    end
  end
  @ode_postamble
end

function ode_solve{uType<:AbstractArray,uEltype,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5}(integrator::ODEIntegrator{Rosenbrock32,uType,uEltype,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  local k₁::uType = similar(u)
  local k₂ = similar(u)
  local k₃::uType = similar(u)
  local tmp::uType
  const kshortsize = 1
  function vecf(t,u,du)
    f(t,reshape(u,sizeu...),reshape(du,sizeu...))
    u = vec(u)
    du = vec(du)
  end
  function vecfreturn(t,u,du)
    f(t,reshape(u,sizeu...),reshape(du,sizeu...))
    return vec(du)
  end
  du1 = zeros(u)
  du2 = zeros(u)
  f₀ = similar(u)
  f₁ = similar(u)
  vectmp3 = similar(vec(u))
  utmp = similar(u); vectmp2 = similar(vec(u))
  dT = similar(u); vectmp = similar(vec(u))
  J = zeros(uEltype,length(u),length(u))
  W = similar(J); tmp2 = similar(u)
  uidx = eachindex(u)
  jidx = eachindex(J)
  f(t,u,fsalfirst)
  if calck
    k = fsalfirst
  end
  cache = (u,du1,du2,uprev,kprev,f₀,f₁,vectmp3,utmp,vectmp2,dT,vectmp,tmp2,k)
  jaccache = (jidx,J,W)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      ForwardDiff.derivative!(dT,(t)->vecfreturn(t,u,du2),t) # Time derivative
      ForwardDiff.jacobian!(J,(du1,u)->vecf(t,u,du1),vec(du1),vec(u))

      W[:] = I-dt*d*J # Can an allocation be cut here?
      @into! vectmp = W\vec(fsalfirst + dt*d*dT)
      k₁ = reshape(vectmp,sizeu...)
      for i in uidx
        utmp[i]=u[i]+dt*k₁[i]/2
      end
      f(t+dt/2,utmp,f₁)
      @into! vectmp2 = W\vec(f₁-k₁)
      tmp = reshape(vectmp2,sizeu...)
      for i in uidx
        k₂[i] = tmp[i] + k₁[i]
      end
      for i in uidx
        tmp[i] = u[i] + dt*k₂[i]
      end
      f(t+dt,tmp,fsallast)
      @into! vectmp3 = W\vec(fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT)
      k₃ = reshape(vectmp3,sizeu...)
      if adaptive
        for i in uidx
          utmp[i] = u[i] + dt*(k₁[i] + 4k₂[i] + k₃[i])/6
          tmp2[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)/(abstol+u[i]*reltol)
        end
        EEst = internalnorm(tmp2)
      else
        for i in uidx
          u[i] = u[i] + dt*(k₁[i] + 4k₂[i] + k₃[i])/6
        end
        f(t,u,fsallast)
      end
      @ode_loopfooter
      recursivecopy!(fsalfirst,fsallast)
    end
  end
  @ode_postamble
end

function ode_solve{uType<:Number,uEltype,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5}(integrator::ODEIntegrator{Rosenbrock32,uType,uEltype,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4,F5})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  local dT::uType
  local J::uType
  #f₀ = fsalfirst
  local k₁::uType
  local f₁::uType
  #f₂ = fsallast
  local k₂::uType
  local k₃::uType
  local tmp::uType
  const kshortsize = 1
  fsalfirst = f(t,u)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      # Time derivative
      dT = ForwardDiff.derivative((t)->f(t,u),t)
      J = ForwardDiff.derivative((u)->f(t,u),u)
      W = I-dt*d*J
      #f₀ = f(t,u)
      if calck
        k = fsalfirst
      end
      k₁ = W\(fsalfirst + dt*d*dT)
      f₁ = f(t+dt/2,u+dt*k₁/2)
      k₂ = W\(f₁-k₁) + k₁
      tmp = u + dt*k₂
      fsallast = f(t+dt,tmp)
      k₃ = W\(fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT)
      if adaptive
        utmp = u + dt*(k₁ + 4k₂ + k₃)/6
        EEst = abs((dt*(k₁ - 2k₂ + k₃)/6)./(abstol+u*reltol))
      else
        u = u + dt*(k₁ + 4k₂ + k₃)/6
        fsallast = f(t,u)
      end
      @ode_loopfooter
      fsalfirst = fsallast
    end
  end
  @ode_postamble
end
