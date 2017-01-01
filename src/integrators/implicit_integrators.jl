type RHS_IE_Scalar{F,uType,tType} <: Function
  f::F
  u_old::uType
  t::tType
  dt::tType
end

function (p::RHS_IE_Scalar)(u,resid)
  resid[1] = u[1] - p.u_old[1] - p.dt*p.f(p.t+p.dt,u)[1]
end

function ode_solve{uType<:Number,algType<:ImplicitEuler,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O})
  @ode_preamble
  local nlres::NLsolve.SolverResults{uEltype}

  uhold::Vector{uType} = Vector{uType}(1)
  u_old::Vector{uType} = Vector{uType}(1)

  rhs = RHS_IE_Scalar(f,u_old,t,dt)

  uhold[1] = u; u_old[1] = u

  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  end
  @inbounds for T in Ts
    while integrator.tdir*t < integrator.tdir*T
      @ode_loopheader
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
      utmp = uhold[1]
      @ode_loopfooter
    end
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

function (p::RHS_IE)(u,resid)
  du = get_du(p.dual_cache, eltype(u))
  p.f(p.t+p.dt,reshape(u,p.sizeu),du)
  for i in p.uidx
    resid[i] = u[i] - p.u_old[i] - p.dt*du[i]
  end
end

function ode_solve{uType<:AbstractArray,algType<:ImplicitEuler,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O})
  @ode_preamble
  local nlres::NLsolve.SolverResults{uEltype}
  uidx = eachindex(u)
  sizeu = size(u) # Change to dynamic by call overloaded type

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,integrator.uprev,integrator.kprev)
  @unpack u_old,dual_cache = cache

  uhold = vec(utmp)

  rhs = RHS_IE(f,u_old,t,dt,sizeu,dual_cache,uidx)

  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  end

  @inbounds for T in Ts
    while integrator.tdir*t < integrator.tdir*T
      @ode_loopheader
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
        f(t+dt,utmp,k)
      end
      @ode_loopfooter
    end
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

function (p::RHS_Trap)(u,resid)
  du1 = get_du(p.dual_cache, eltype(u)); du2 = get_du(p.dual_cache2, eltype(p.u_old))
  p.f(p.t,reshape(p.u_old,p.sizeu),du2)
  p.f(p.t+p.dt,reshape(u,p.sizeu),du1)
  for i in p.uidx
    resid[i] = u[i] - p.u_old[i] - (p.dt/2)*(du1[i]+du2[i])
  end
end

function ode_solve{uType<:AbstractArray,algType<:Trapezoid,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O})
  @ode_preamble
  local nlres::NLsolve.SolverResults{uEltype}
  uidx = eachindex(u)
  sizeu = size(u) # Change to dynamic by call overloaded type

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,integrator.uprev,integrator.kprev)
  @unpack u_old,dual_cache,dual_cache2 = cache

  dto2 = dt/2

  rhs = RHS_Trap(f,u_old,t,dt,sizeu,dual_cache,dual_cache2,uidx)

  uhold = vec(utmp)

  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  end

  @inbounds for T in Ts
    while integrator.tdir*t < integrator.tdir*T
      @ode_loopheader
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
        f(t+dt,utmp,k)
      end
      @ode_loopfooter
    end
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

function (p::RHS_Trap_Scalar)(u,resid)
  resid[1] = u[1] - p.u_old[1] - (p.dt/2)*(p.f(p.t,p.u_old)[1] + p.f(p.t+p.dt,u)[1])
end

function ode_solve{uType<:Number,algType<:Trapezoid,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O}(integrator::ODEIntegrator{algType,uType,tType,tTypeNoUnits,ksEltype,SolType,rateType,F,ECType,O})
  @ode_preamble
  dto2::tType = dt/2
  local nlres::NLsolve.SolverResults{uEltype}
  uhold::Vector{uType} = Vector{uType}(1)
  u_old::Vector{uType} = Vector{uType}(1)
  uhold[1] = u; u_old[1] = u

  rhs = RHS_Trap_Scalar(f,u_old,t,dt)

  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  end

  @inbounds for T in Ts
      while integrator.tdir*t < integrator.tdir*T
      @ode_loopheader
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
      utmp = uhold[1]
      @ode_loopfooter
    end
  end
  ode_postamble!(integrator)
  nothing
end
