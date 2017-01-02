function ode_solve{uType<:AbstractArray,algType<:Rosenbrock23,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))


  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp2 = integrator.cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast
  sizeu = size(uprev) # Change to dynamic by call overloaded type
  integrator.kshortsize = 2
  function vecf(t,uprev,du)
    f(t,reshape(uprev,sizeu...),reshape(du,sizeu...))
    uprev = vec(uprev)
    du = vec(du)
  end
  function vecfreturn(t,uprev,du)
    f(t,reshape(uprev,sizeu...),reshape(du,sizeu...))
    return vec(du)
  end
  tmp = reshape(vectmp2,sizeu...)
  uidx = eachindex(uprev)
  jidx = eachindex(J)
  integrator.k = [k₁,k₂]
  f(t,uprev,fsalfirst)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      #if alg_autodiff(alg)
        ForwardDiff.derivative!(dT,(t)->vecfreturn(t,uprev,du2),t) # Time derivative of each component
        ForwardDiff.jacobian!(J,(du1,uprev)->vecf(t,uprev,du1),vec(du1),vec(uprev))
      #else
      #  Calculus.finite_difference!((t)->vecfreturn(t,uprev,du2),[t],dT)
      #  Calculus.finite_difference_jacobian!((du1,uprev)->vecf(t,uprev,du1),vec(uprev),vec(du1),J)
      #end

      W[:] = I-dt*d*J # Can an allocation be cut here?
      @into! vectmp = W\vec(fsalfirst + dt*d*dT)
      recursivecopy!(k₁,reshape(vectmp,sizeu...))
      for i in uidx
        u[i]=uprev[i]+dt*k₁[i]/2
      end
      f(t+dt/2,u,f₁)
      @into! vectmp2 = W\vec(f₁-k₁)
      for i in uidx
        k₂[i] = tmp[i] + k₁[i]
        u[i] = uprev[i] + dt*k₂[i]
      end
      if integrator.opts.adaptive
        f(t+dt,u,integrator.fsallast)
        @into! vectmp3 = W\vec(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+dt*dT)
        k₃ = reshape(vectmp3,sizeu...)
        for i in uidx
          tmp2[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)./(integrator.opts.abstol+max(abs(uprev[i]),abs(u[i]))*integrator.opts.reltol)
        end
        integrator.EEst = integrator.opts.internalnorm(tmp2)
      end
      @pack_integrator
      ode_loopfooter!(integrator)
      @unpack_integrator
      if isempty(integrator.tstops)
        break
      end
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,algType<:Rosenbrock23,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
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
  integrator.k = k
  integrator.fsalfirst = f(t,uprev)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      # Time derivative
      dT = ForwardDiff.derivative((t)->f(t,uprev),t)
      J = ForwardDiff.derivative((uprev)->f(t,uprev),uprev)
      W = 1-dt*d*J
      k₁ = W\(integrator.fsalfirst + dt*d*dT)
      f₁ = f(t+dt/2,uprev+dt*k₁/2)
      k₂ = W\(f₁-k₁) + k₁
      u = uprev + dt*k₂
      if integrator.opts.adaptive
        integrator.fsallast = f(t+dt,u)
        k₃ = W\(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT)
        integrator.EEst = abs((dt*(k₁ - 2k₂ + k₃)/6)./(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
      end
      if integrator.opts.calck
        k[1] = k₁
        k[2] = k₂
      end
      @pack_integrator
      ode_loopfooter!(integrator)
      @unpack_integrator
      if isempty(integrator.tstops)
        break
      end
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,algType<:Rosenbrock32,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))


  @unpack k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp2 = integrator.cache

  integrator.fsalfirst = fsalfirst
  integrator.fsallast = fsallast

  sizeu = size(uprev) # Change to dynamic by call overloaded type
  integrator.kshortsize = 2
  function vecf(t,uprev,du)
    f(t,reshape(uprev,sizeu...),reshape(du,sizeu...))
    uprev = vec(uprev)
    du = vec(du)
  end
  function vecfreturn(t,uprev,du)
    f(t,reshape(uprev,sizeu...),reshape(du,sizeu...))
    return vec(du)
  end
  uidx = eachindex(uprev)
  jidx = eachindex(J)
  integrator.k = [k₁,k₂]
  f(t,uprev,integrator.fsalfirst)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      ForwardDiff.derivative!(dT,(t)->vecfreturn(t,uprev,du2),t) # Time derivative
      ForwardDiff.jacobian!(J,(du1,uprev)->vecf(t,uprev,du1),vec(du1),vec(uprev))

      W[:] = I-dt*d*J # Can an allocation be cut here?
      @into! vectmp = W\vec(integrator.fsalfirst + dt*d*dT)
      recursivecopy!(k₁,reshape(vectmp,sizeu...))
      for i in uidx
        u[i]=uprev[i]+dt*k₁[i]/2
      end
      f(t+dt/2,u,f₁)
      @into! vectmp2 = W\vec(f₁-k₁)
      tmp = reshape(vectmp2,sizeu...)
      for i in uidx
        k₂[i] = tmp[i] + k₁[i]
      end
      for i in uidx
        tmp[i] = uprev[i] + dt*k₂[i]
      end
      f(t+dt,tmp,integrator.fsallast)
      @into! vectmp3 = W\vec(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT)
      k₃ = reshape(vectmp3,sizeu...)
      for i in uidx
        u[i] = uprev[i] + dt*(k₁[i] + 4k₂[i] + k₃[i])/6
      end
      if integrator.opts.adaptive
        for i in uidx
          tmp2[i] = (dt*(k₁[i] - 2k₂[i] + k₃[i])/6)/(integrator.opts.abstol+max(abs(uprev[i]),abs(u[i]))*integrator.opts.reltol)
        end
        integrator.EEst = integrator.opts.internalnorm(tmp2)
      end
      @pack_integrator
      ode_loopfooter!(integrator)
      @unpack_integrator
      if isempty(integrator.tstops)
        break
      end
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,algType<:Rosenbrock32,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
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
  k = ksEltype(2)
  integrator.k = k
  integrator.fsalfirst = f(t,uprev)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      # Time derivative
      dT = ForwardDiff.derivative((t)->f(t,uprev),t)
      J = ForwardDiff.derivative((uprev)->f(t,uprev),uprev)
      W = 1-dt*d*J
      #f₀ = f(t,uprev)
      k₁ = W\(integrator.fsalfirst + dt*d*dT)
      f₁ = f(t+dt/2,uprev+dt*k₁/2)
      k₂ = W\(f₁-k₁) + k₁
      tmp = uprev + dt*k₂
      integrator.fsallast = f(t+dt,tmp)
      k₃ = W\(integrator.fsallast - c₃₂*(k₂-f₁)-2(k₁-integrator.fsalfirst)+dt*dT)
      u = uprev + dt*(k₁ + 4k₂ + k₃)/6
      if integrator.opts.adaptive
        integrator.EEst = abs((dt*(k₁ - 2k₂ + k₃)/6)./(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
      end
      if integrator.opts.calck
        k[1] = k₁
        k[2] = k₂
      end
      @pack_integrator
      ode_loopfooter!(integrator)
      @unpack_integrator
      if isempty(integrator.tstops)
        break
      end
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end
