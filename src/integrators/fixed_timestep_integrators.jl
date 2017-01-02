function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Euler,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  integrator.fsalfirst = f(t,u) # For the interpolation, needs k at the updated point
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      k = integrator.fsalfirst
      utmp = muladd(dt,k,u)
      k = f(t+dt,utmp) # For the interpolation, needs k at the updated point
      integrator.fsallast = k
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Euler,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  uidx = eachindex(u)

  @unpack k,fsalfirst = integrator.cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = k
  f(t,u,integrator.fsalfirst) # For the interpolation, needs k at the updated point
  @inbounds while !isempty(integrator.tstops)
      while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      for i in uidx
        utmp[i] = muladd(dt,integrator.fsalfirst[i],u[i])
      end
      f(t+dt,utmp,integrator.fsallast) # For the interpolation, needs k at the updated point
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Midpoint,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  halfdt::tType = dt/2
  local du::rateType
  integrator.fsalfirst = f(t,u)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      k = integrator.fsalfirst
      k = f(t+halfdt,u+halfdt*k)
      utmp = u + dt*k
      integrator.fsallast = f(t+dt,utmp) # For interpolation, then FSAL'd
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Midpoint,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  halfdt::tType = dt/2
  uidx = eachindex(u)
  if integrator.opts.calck # Not initialized if not dense
    if integrator.calcprevs
      integrator.kprev = similar(rate_prototype)
    end
  end


  @unpack k,du,utilde,fsalfirst = integrator.cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = k
  f(t,u,integrator.fsalfirst) # FSAL for interpolation
  @inbounds while !isempty(integrator.tstops)
      while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      for i in uidx
        utilde[i] = muladd(halfdt,integrator.fsalfirst[i],u[i])
      end
      f(t+halfdt,utilde,du)
      for i in uidx
        utmp[i] = muladd(dt,du[i],u[i])
      end
      f(t+dt,utmp,k)
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end

function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{RK4,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  halfdt::tType = dt/2
  local k₁::rateType
  local k₂::rateType
  local k₃::rateType
  local k₄::rateType
  local ttmp::tType
  integrator.fsalfirst = f(t,u)
  @inbounds while !isempty(integrator.tstops)
      while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      k₁=integrator.fsalfirst
      ttmp = t+halfdt
      k₂ = f(ttmp,muladd(halfdt,k₁,u))
      k₃ = f(ttmp,muladd(halfdt,k₂,u))
      k₄ = f(t+dt,muladd(dt,k₃,u))
      utmp = muladd(dt/6,muladd(2,(k₂ + k₃),k₁+k₄),u)
      k = f(t+dt,utmp)
      integrator.fsallast = k
      @ode_loopfooter
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

  uidx = eachindex(u)


  @unpack tmp,k₁,k₂,k₃,k₄,k = integrator.cache
  integrator.fsalfirst = k₁
  integrator.fsallast = k
  integrator.k = k
  f(t,u,k₁) # pre-start FSAL
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*t < integrator.tdir*top(integrator.tstops)
      @ode_loopheader
      ttmp = t+halfdt
      for i in uidx
        tmp[i] = muladd(halfdt,k₁[i],u[i])
      end
      f(ttmp,tmp,k₂)
      for i in uidx
        tmp[i] = muladd(halfdt,k₂[i],u[i])
      end
      f(ttmp,tmp,k₃)
      for i in uidx
        tmp[i] = muladd(dt,k₃[i],u[i])
      end
      f(t+dt,tmp,k₄)
      for i in uidx
        utmp[i] = muladd(dt/6,muladd(2,(k₂[i] + k₃[i]),k₁[i] + k₄[i]),u[i])
      end
      f(t+dt,utmp,k)
      @ode_loopfooter
    end
    !isempty(integrator.tstops) && pop!(integrator.tstops)
  end
  ode_postamble!(integrator)
  nothing
end
