function ode_solve{uType<:AbstractArray,algType<:Rosenbrock23,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,integrator.uprev,integrator.kprev)
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp2 = cache

  sizeu = size(u) # Change to dynamic by call overloaded type
  integrator.kshortsize = 2
  function vecf(t,u,du)
    f(t,reshape(u,sizeu...),reshape(du,sizeu...))
    u = vec(u)
    du = vec(du)
  end
  function vecfreturn(t,u,du)
    f(t,reshape(u,sizeu...),reshape(du,sizeu...))
    return vec(du)
  end
  tmp = reshape(vectmp2,sizeu...)
  uidx = eachindex(u)
  jidx = eachindex(J)
  k = [k₁,k₂]
  f(t,u,fsalfirst)
  @inbounds for T in Ts
    while integrator.tdir*t < integrator.tdir*T
      @ode_loopheader
      #if alg_autodiff(alg)
        ForwardDiff.derivative!(dT,(t)->vecfreturn(t,u,du2),t) # Time derivative of each component
        ForwardDiff.jacobian!(J,(du1,u)->vecf(t,u,du1),vec(du1),vec(u))
      #else
      #  Calculus.finite_difference!((t)->vecfreturn(t,u,du2),[t],dT)
      #  Calculus.finite_difference_jacobian!((du1,u)->vecf(t,u,du1),vec(u),vec(du1),J)
      #end

      W[:] = I-dt*d*J # Can an allocation be cut here?
      @into! vectmp = W\vec(fsalfirst + dt*d*dT)
      recursivecopy!(k₁,reshape(vectmp,sizeu...))
      for i in uidx
        utmp[i]=u[i]+dt*k₁[i]/2
      end
      f(t+dt/2,utmp,f₁)
      @into! vectmp2 = W\vec(f₁-k₁)
      for i in uidx
        k₂[i] = tmp[i] + k₁[i]
        utmp[i] = u[i] + dt*k₂[i]
      end
      if integrator.opts.adaptive
        f(t+dt,utmp,fsallast)
        @into! vectmp3 = W\vec(fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT)
        k₃ = reshape(vectmp3,sizeu...)
        for i in uidx
          tmp2[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)./(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol)
        end
        EEst = integrator.opts.internalnorm(tmp2)
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,algType<:Rosenbrock23,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  local dT::uType
  local J::uType
  local k₁::uType
  local f₁::uType
  local k₂::uType
  local k₃::uType
  integrator.kshortsize = 2
  k = ksEltype(2)
  fsalfirst = f(t,u)
  @inbounds for T in Ts
    while integrator.tdir*t < integrator.tdir*T
      @ode_loopheader
      # Time derivative
      dT = ForwardDiff.derivative((t)->f(t,u),t)
      J = ForwardDiff.derivative((u)->f(t,u),u)
      W = 1-dt*d*J
      k₁ = W\(fsalfirst + dt*d*dT)
      f₁ = f(t+dt/2,u+dt*k₁/2)
      k₂ = W\(f₁-k₁) + k₁
      utmp = u + dt*k₂
      if integrator.opts.adaptive
        fsallast = f(t+dt,utmp)
        k₃ = W\(fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT)
        EEst = abs((dt*(k₁ - 2k₂ + k₃)/6)./(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol))
      end
      if integrator.opts.calck
        k[1] = k₁
        k[2] = k₂
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,algType<:Rosenbrock32,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,integrator.uprev,integrator.kprev)
  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp2 = cache

  sizeu = size(u) # Change to dynamic by call overloaded type
  integrator.kshortsize = 2
  function vecf(t,u,du)
    f(t,reshape(u,sizeu...),reshape(du,sizeu...))
    u = vec(u)
    du = vec(du)
  end
  function vecfreturn(t,u,du)
    f(t,reshape(u,sizeu...),reshape(du,sizeu...))
    return vec(du)
  end
  uidx = eachindex(u)
  jidx = eachindex(J)
  k = [k₁,k₂]
  f(t,u,fsalfirst)
  @inbounds for T in Ts
    while integrator.tdir*t < integrator.tdir*T
      @ode_loopheader
      ForwardDiff.derivative!(dT,(t)->vecfreturn(t,u,du2),t) # Time derivative
      ForwardDiff.jacobian!(J,(du1,u)->vecf(t,u,du1),vec(du1),vec(u))

      W[:] = I-dt*d*J # Can an allocation be cut here?
      @into! vectmp = W\vec(fsalfirst + dt*d*dT)
      recursivecopy!(k₁,reshape(vectmp,sizeu...))
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
      for i in uidx
        utmp[i] = u[i] + dt*(k₁[i] + 4k₂[i] + k₃[i])/6
      end
      if integrator.opts.adaptive
        for i in uidx
          tmp2[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)/(integrator.opts.abstol+max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol)
        end
        EEst = integrator.opts.internalnorm(tmp2)
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,algType<:Rosenbrock32,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O})
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
  integrator.kshortsize = 2
  if integrator.opts.calck
    k = ksEltype()
    for i in 1:2
      push!(k,zero(rateType))
    end
    if integrator.calcprevs
      integrator.kprev = deepcopy(k)
      for i in 1:2
        push!(integrator.kprev,zero(rateType))
      end
    end
  end
  fsalfirst = f(t,u)
  @inbounds for T in Ts
    while integrator.tdir*t < integrator.tdir*T
      @ode_loopheader
      # Time derivative
      dT = ForwardDiff.derivative((t)->f(t,u),t)
      J = ForwardDiff.derivative((u)->f(t,u),u)
      W = 1-dt*d*J
      #f₀ = f(t,u)
      k₁ = W\(fsalfirst + dt*d*dT)
      f₁ = f(t+dt/2,u+dt*k₁/2)
      k₂ = W\(f₁-k₁) + k₁
      tmp = u + dt*k₂
      fsallast = f(t+dt,tmp)
      k₃ = W\(fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT)
      utmp = u + dt*(k₁ + 4k₂ + k₃)/6
      if integrator.opts.adaptive
        EEst = abs((dt*(k₁ - 2k₂ + k₃)/6)./(integrator.opts.abstol+max(abs(u),abs(utmp))*integrator.opts.reltol))
      end
      if integrator.opts.calck
        k[1] = k₁
        k[2] = k₂
      end
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end
