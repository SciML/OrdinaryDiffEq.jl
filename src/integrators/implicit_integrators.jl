type RHS_IE_Scalar{F,uType,tType} <: Function
  f::F
  u_old::uType
  t::tType
  dt::tType
end

function (p::RHS_IE_Scalar)(u,resid)
  resid[1] = u[1] - p.u_old[1] - p.dt*p.f(p.t+p.dt,u)[1]
end

@inline function initialize!(integrator,cache::ImplicitEulerConstantCache)
  cache.uhold[1] = integrator.uprev; cache.u_old[1] = integrator.uprev
end


function ode_solve{uType<:Number,algType<:ImplicitEuler,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack t,dt,uprev,u,f,k = integrator
      @unpack uhold,u_old,rhs,adf = integrator.cache
      u_old[1] = uhold[1]
      rhs.t = t
      rhs.dt = dt
      if alg_autodiff(integrator.alg)
        nlres = NLsolve.nlsolve(adf,uhold)
      else
        nlres = NLsolve.nlsolve(rhs,uhold,autodiff=alg_autodiff(integrator.alg))
      end
      uhold[1] = nlres.zero[1]
      if integrator.opts.calck
        k = f(t+dt,uhold[1])
      end
      u = uhold[1]
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

type RHS_IE{F,uType,tType,DiffCacheType,SizeType,uidxType} <: Function
  f::F
  u_old::uType
  t::tType
  dt::tType
  dual_cache::DiffCacheType
  sizeu::SizeType
  uidx::uidxType
end
function (p::RHS_IE)(uprev,resid)
  du = get_du(p.dual_cache, eltype(uprev))
  p.f(p.t+p.dt,reshape(uprev,p.sizeu),du)
  for i in p.uidx
    resid[i] = uprev[i] - p.u_old[i] - p.dt*du[i]
  end
end

@inline function initialize!(integrator,cache::ImplicitEulerCache)
  integrator.k = cache.k
end

function ode_solve{uType<:AbstractArray,algType<:ImplicitEuler,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack t,dt,uprev,u,f,k = integrator
      uidx = eachindex(integrator.uprev)
      @unpack u_old,dual_cache,k,adf,rhs,uhold = integrator.cache
      copy!(u_old,uhold)
      rhs.t = t
      rhs.dt = dt
      rhs.uidx = uidx
      rhs.sizeu = size(u)
      if alg_autodiff(integrator.alg)
        nlres = NLsolve.nlsolve(adf,uhold)
      else
        nlres = NLsolve.nlsolve(rhs,uhold,autodiff=alg_autodiff(integrator.alg))
      end
      copy!(uhold,nlres.zero)
      if integrator.opts.calck
        f(t+dt,u,k)
      end
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


type RHS_Trap{F,uType,rateType,tType,SizeType,DiffCacheType,uidxType} <: Function
  f::F
  u_old::uType
  f_old::rateType
  t::tType
  dt::tType
  sizeu::SizeType
  dual_cache::DiffCacheType
  uidx::uidxType
end

function (p::RHS_Trap)(uprev,resid)
  du1 = get_du(p.dual_cache, eltype(uprev))
  p.f(p.t+p.dt,reshape(uprev,p.sizeu),du1)
  for i in p.uidx
    resid[i] = uprev[i] - p.u_old[i] - (p.dt/2)*(du1[i]+p.f_old[i])
  end
end

@inline function initialize!(integrator,cache::TrapezoidCache)
  integrator.k = cache.k
  @unpack k,f_old = integrator.cache
  integrator.fsalfirst = f_old
  integrator.fsallast = cache.k
  integrator.k = k
  integrator.f(integrator.t,integrator.uprev,integrator.fsalfirst)
end

function ode_solve{uType<:AbstractArray,algType<:Trapezoid,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack t,dt,uprev,u,f,k = integrator
      uidx = eachindex(integrator.uprev)
      @unpack u_old,dual_cache,k,rhs,adf,uhold = integrator.cache
      copy!(u_old,uhold)
      # copy!(rhs.f_old,f_old) Implicitly done by pointers: fsalfirst === f_old == rhs.f_old
      rhs.t = t
      rhs.dt = dt
      rhs.uidx = uidx
      rhs.sizeu = size(u)
      if alg_autodiff(integrator.alg)
        nlres = NLsolve.nlsolve(adf,uhold)
      else
        nlres = NLsolve.nlsolve(rhs,uhold,autodiff=alg_autodiff(integrator.alg))
      end
      copy!(uhold,nlres.zero)
      if integrator.opts.calck
        f(t+dt,u,k)
      end
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

type RHS_Trap_Scalar{F,uType,rateType,tType} <: Function
  f::F
  u_old::uType
  f_old::rateType
  t::tType
  dt::tType
end

function (p::RHS_Trap_Scalar)(uprev,resid)
  resid[1] = uprev[1] - p.u_old[1] - (p.dt/2)*(p.f_old + p.f(p.t+p.dt,uprev)[1])
end

@inline function initialize!(integrator,cache::TrapezoidConstantCache)
  cache.uhold[1] = integrator.uprev; cache.u_old[1] = integrator.uprev
  integrator.fsalfirst = integrator.f(integrator.t,integrator.uprev)
end

function ode_solve{uType<:Number,algType<:Trapezoid,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.tstops)
      while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack t,dt,uprev,u,f,k = integrator
      @unpack uhold,u_old,rhs,adf = integrator.cache
      u_old[1] = uhold[1]
      rhs.t = t
      rhs.dt = dt
      rhs.f_old = integrator.fsalfirst
      if alg_autodiff(integrator.alg)
        nlres = NLsolve.nlsolve(adf,uhold)
      else
        nlres = NLsolve.nlsolve(rhs,uhold,autodiff=alg_autodiff(integrator.alg))
      end
      uhold[1] = nlres.zero[1]
      k = f(t+dt,uhold[1])
      integrator.fsallast = k
      u = uhold[1]
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
