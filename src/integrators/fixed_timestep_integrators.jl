@inline function initialize!(integrator,cache::EulerConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
end

function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Euler,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack t,dt,uprev,u,f,k = integrator
      k = integrator.fsalfirst
      u = muladd(dt,k,uprev)
      k = f(t+dt,u) # For the interpolation, needs k at the updated point
      integrator.fsallast = k
      @pack integrator = t,dt,u,k
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

@inline function initialize!(integrator,cache::EulerCache)
  @unpack k,fsalfirst = integrator.cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = k
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Euler,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
      while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack t,dt,uprev,u,f,k = integrator
      uidx = eachindex(integrator.uprev)
      for i in uidx
        u[i] = muladd(dt,integrator.fsalfirst[i],uprev[i])
      end
      f(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
      @pack integrator = t,dt,u,k
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

@inline function initialize!(integrator,cache::MidpointConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
end

function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Midpoint,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack t,dt,uprev,u,f,k = integrator
            halfdt = dt/2
      k = integrator.fsalfirst
      k = f(t+halfdt,uprev+halfdt*k)
      u = uprev + dt*k
      integrator.fsallast = f(t+dt,u) # For interpolation, then FSAL'd
      @pack integrator = t,dt,u,k
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

@inline function initialize!(integrator,cache::MidpointCache)
  if integrator.opts.calck # Not initialized if not dense
    if integrator.calcprevs
      integrator.kprev = similar(integrator.rate_prototype)
    end
  end
  @unpack k,fsalfirst = integrator.cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = k
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # FSAL for interpolation
end

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{Midpoint,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
      while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack t,dt,uprev,u,f,k = integrator
      uidx = eachindex(integrator.uprev)
      @unpack k,du,utilde,fsalfirst = integrator.cache
      halfdt::tType = dt/2
      for i in uidx
        utilde[i] = muladd(halfdt,integrator.fsalfirst[i],uprev[i])
      end
      f(t+halfdt,utilde,du)
      for i in uidx
        u[i] = muladd(dt,du[i],uprev[i])
      end
      f(t+dt,u,k)
      @pack integrator = t,dt,u,k
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

@inline function initialize!(integrator,cache::RK4ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev) # Pre-start fsal
end

function ode_solve{uType<:Number,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{RK4,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
      while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack t,dt,uprev,u,f,k = integrator
            halfdt = dt/2
      k₁ =integrator.fsalfirst
      ttmp = t+halfdt
      k₂ = f(ttmp,muladd(halfdt,k₁,uprev))
      k₃ = f(ttmp,muladd(halfdt,k₂,uprev))
      k₄ = f(t+dt,muladd(dt,k₃,uprev))
      u = muladd(dt/6,muladd(2,(k₂ + k₃),k₁+k₄),uprev)
      k = f(t+dt,u)
      integrator.fsallast = k
      @pack integrator = t,dt,u,k
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

@inline function initialize!(integrator,cache::RK4Cache)
  if integrator.calcprevs
    integrator.kprev = similar(integrator.rate_prototype)
  end
  @unpack tmp,k₁,k₂,k₃,k₄,k = integrator.cache
  integrator.fsalfirst = k₁
  integrator.fsallast = k
  integrator.k = k
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst) # pre-start FSAL
end

function ode_solve{uType<:AbstractArray,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{RK4,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack t,dt,uprev,u,f,k = integrator
      uidx = eachindex(integrator.uprev)
      @unpack tmp,k₁,k₂,k₃,k₄,k = integrator.cache
      halfdt = dt/2
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
      @pack integrator = t,dt,u,k
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
