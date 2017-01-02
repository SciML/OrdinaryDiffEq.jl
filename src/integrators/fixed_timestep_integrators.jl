function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Euler,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  integrator.fsalfirst = f(t,uprev) # For the interpolation, needs k at the updated point
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack_integrator
      k = integrator.fsalfirst
      u = muladd(dt,k,uprev)
      k = f(t+dt,u) # For the interpolation, needs k at the updated point
      integrator.fsallast = k
      @pack_integrator
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

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Euler,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  uidx = eachindex(uprev)

  @unpack k,fsalfirst = integrator.cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = k
  f(t,uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
  @inbounds while !isempty(integrator.tstops)
      while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack_integrator
      for i in uidx
        u[i] = muladd(dt,integrator.fsalfirst[i],uprev[i])
      end
      f(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
      @pack_integrator
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

function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Midpoint,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  integrator.fsalfirst = f(t,uprev)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack_integrator
      halfdt = dt/2
      k = integrator.fsalfirst
      k = f(t+halfdt,uprev+halfdt*k)
      u = uprev + dt*k
      integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
      @pack_integrator
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

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Midpoint,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  halfdt::tType = dt/2
  uidx = eachindex(uprev)
  if integrator.opts.calck # Not initialized if not dense
    if integrator.calcprevs
      integrator.kprev = similar(rate_prototype)
    end
  end


  @unpack k,du,utilde,fsalfirst = integrator.cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = k
  f(t,uprev,integrator.fsalfirst) # FSAL for interpolation
  @inbounds while !isempty(integrator.tstops)
      while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack_integrator
      for i in uidx
        utilde[i] = muladd(halfdt,integrator.fsalfirst[i],uprev[i])
      end
      f(t+halfdt,utilde,du)
      for i in uidx
        u[i] = muladd(dt,du[i],uprev[i])
      end
      f(t+dt,u,k)
      @pack_integrator
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

function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{RK4,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  integrator.fsalfirst = f(t,uprev)
  @inbounds while !isempty(integrator.tstops)
      while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack_integrator
      halfdt = dt/2
      k₁=integrator.fsalfirst
      ttmp = t+halfdt
      k₂ = f(ttmp,muladd(halfdt,k₁,uprev))
      k₃ = f(ttmp,muladd(halfdt,k₂,uprev))
      k₄ = f(t+dt,muladd(dt,k₃,uprev))
      u = muladd(dt/6,muladd(2,(k₂ + k₃),k₁+k₄),uprev)
      k = f(t+dt,u)
      integrator.fsallast = k
      @pack_integrator
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

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{RK4,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  halfdt::tType = dt/2

  if integrator.calcprevs
    integrator.kprev = similar(rate_prototype)
  end

  uidx = eachindex(uprev)


  @unpack tmp,k₁,k₂,k₃,k₄,k = integrator.cache
  integrator.fsalfirst = k₁
  integrator.fsallast = k
  integrator.k = k
  f(t,uprev,k₁) # pre-start FSAL
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack_integrator
      ttmp = t+halfdt
      for i in uidx
        tmp[i] = muladd(halfdt,k₁[i],uprev[i])
      end
      f(ttmp,tmp,k₂)
      for i in uidx
        tmp[i] = muladd(halfdt,k₂[i],uprev[i])
      end
      f(ttmp,tmp,k₃)
      for i in uidx
        tmp[i] = muladd(dt,k₃[i],uprev[i])
      end
      f(t+dt,tmp,k₄)
      for i in uidx
        u[i] = muladd(dt/6,muladd(2,(k₂[i] + k₃[i]),k₁[i] + k₄[i]),uprev[i])
      end
      f(t+dt,u,k)
      @pack_integrator
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
