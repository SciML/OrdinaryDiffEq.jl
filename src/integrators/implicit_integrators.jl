type RHS_IE_Scalar{F,uType,tType} <: Function
  f::F
  u_old::uType
  t::tType
  dt::tType
end

function (p::RHS_IE_Scalar)(uprev,resid)
  resid[1] = uprev[1] - p.u_old[1] - p.dt*p.f(p.t+p.dt,uprev)[1]
end

function ode_solve{uType<:Number,algType<:ImplicitEuler,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  uhold = Vector{uType}(1)
  u_old = Vector{uType}(1)

  rhs = RHS_IE_Scalar(f,u_old,t,dt)

  uhold[1] = uprev; u_old[1] = uprev

  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,integrator.cache)
  end
  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack_integrator
      u_old[1] = uhold[1]
      rhs.t = t
      rhs.dt = dt
      if alg_autodiff(alg)
        nlres = NLsolve.nlsolve(adf,uhold)
      else
        nlres = NLsolve.nlsolve(rhs,uhold,autodiff=alg_autodiff(alg))
      end
      uhold[1] = nlres.zero[1]
      if integrator.opts.calck
        k = f(t+dt,uhold[1])
      end
      u = uhold[1]
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

type RHS_IE{F,uType,tType,SizeType,DiffCacheType,uidxType} <: Function
  f::F
  u_old::uType
  t::tType
  dt::tType
  sizeu::SizeType
  dual_cache::DiffCacheType
  uidx::uidxType
end

function (p::RHS_IE)(uprev,resid)
  du = get_du(p.dual_cache, eltype(uprev))
  p.f(p.t+p.dt,reshape(uprev,p.sizeu),du)
  for i in p.uidx
    resid[i] = uprev[i] - p.u_old[i] - p.dt*du[i]
  end
end

function ode_solve{uType<:AbstractArray,algType<:ImplicitEuler,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  uidx = eachindex(uprev)
  sizeu = size(uprev) # Change to dynamic by call overloaded type


  @unpack u_old,dual_cache,k = integrator.cache
  integrator.k = k

  uhold = vec(u)

  rhs = RHS_IE(f,u_old,t,dt,sizeu,dual_cache,uidx)

  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,integrator.cache)
  end

  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack_integrator
      copy!(u_old,uhold)
      rhs.t = t
      rhs.dt = dt
      # rhs.uidx = uidx
      # rhs.sizeu = sizeu
      if alg_autodiff(alg)
        nlres = NLsolve.nlsolve(adf,uhold)
      else
        nlres = NLsolve.nlsolve(rhs,uhold,autodiff=alg_autodiff(alg))
      end
      copy!(uhold,nlres.zero)
      if integrator.opts.calck
        f(t+dt,u,k)
      end
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


type RHS_Trap{F,uType,tType,SizeType,DiffCacheType,uidxType} <: Function
  f::F
  u_old::uType
  t::tType
  dt::tType
  sizeu::SizeType
  dual_cache::DiffCacheType
  dual_cache2::DiffCacheType
  uidx::uidxType
end

function (p::RHS_Trap)(uprev,resid)
  du1 = get_du(p.dual_cache, eltype(uprev)); du2 = get_du(p.dual_cache2, eltype(p.u_old))
  p.f(p.t,reshape(p.u_old,p.sizeu),du2)
  p.f(p.t+p.dt,reshape(uprev,p.sizeu),du1)
  for i in p.uidx
    resid[i] = uprev[i] - p.u_old[i] - (p.dt/2)*(du1[i]+du2[i])
  end
end

function ode_solve{uType<:AbstractArray,algType<:Trapezoid,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  uidx = eachindex(uprev)
  sizeu = size(uprev) # Change to dynamic by call overloaded type


  @unpack u_old,dual_cache,dual_cache2,k = integrator.cache
  integrator.k = k

  dto2 = dt/2

  rhs = RHS_Trap(f,u_old,t,dt,sizeu,dual_cache,dual_cache2,uidx)

  uhold = vec(u)

  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,integrator.cache)
  end

  @inbounds while !isempty(integrator.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack_integrator
      copy!(u_old,uhold)
      rhs.t = t
      rhs.dt = dt
      # rhs.uidx = uidx
      # rhs.sizeu = sizeu
      if alg_autodiff(alg)
        nlres = NLsolve.nlsolve(adf,uhold)
      else
        nlres = NLsolve.nlsolve(rhs,uhold,autodiff=alg_autodiff(alg))
      end
      copy!(uhold,nlres.zero)
      if integrator.opts.calck
        f(t+dt,u,k)
      end
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

type RHS_Trap_Scalar{F,uType,tType} <: Function
  f::F
  u_old::uType
  t::tType
  dt::tType
end

function (p::RHS_Trap_Scalar)(uprev,resid)
  resid[1] = uprev[1] - p.u_old[1] - (p.dt/2)*(p.f(p.t,p.u_old)[1] + p.f(p.t+p.dt,uprev)[1])
end

function ode_solve{uType<:Number,algType<:Trapezoid,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O})
  @ode_preamble
  dto2::tType = dt/2
  uhold::Vector{uType} = Vector{uType}(1)
  u_old::Vector{uType} = Vector{uType}(1)
  uhold[1] = uprev; u_old[1] = uprev

  rhs = RHS_Trap_Scalar(f,u_old,t,dt)

  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,integrator.cache)
  end

  @inbounds while !isempty(integrator.tstops)
      while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      @ode_exit_conditions
      @unpack_integrator
      u_old[1] = uhold[1]
      rhs.t = t
      rhs.dt = dt
      if alg_autodiff(alg)
        nlres = NLsolve.nlsolve(adf,uhold)
      else
        nlres = NLsolve.nlsolve(rhs,uhold,autodiff=alg_autodiff(alg))
      end
      uhold[1] = nlres.zero[1]
      if integrator.opts.calck
        k = f(t+dt,uhold[1])
      end
      u = uhold[1]
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
