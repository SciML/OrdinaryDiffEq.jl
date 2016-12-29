type DEOptions{uEltype,uEltypeNoUnits,tTypeNoUnits,tType,F2,F3,F4,F5}
  maxiters::Int
  timeseries_steps::Int
  save_timeseries::Bool
  adaptive::Bool
  abstol::uEltype
  reltol::uEltypeNoUnits
  gamma::tTypeNoUnits
  qmax::tTypeNoUnits
  qmin::tTypeNoUnits
  dtmax::tType
  dtmin::tType
  internalnorm::F2
  progress::Bool
  progress_steps::Int
  progress_name::String
  progress_message::F5
  beta1::tTypeNoUnits
  beta2::tTypeNoUnits
  qoldinit::tTypeNoUnits
  dense::Bool
  saveat::Vector{tType}
  callback::F3
  isoutofdomain::F4
  calck::Bool
end

type ODEIntegrator{algType<:OrdinaryDiffEqAlgorithm,uType<:Union{AbstractArray,Number},tType,ksEltype,SolType,rateType,F,ECType,O}
  sol::SolType
  u::uType
  k::ksEltype
  t::tType
  dt::tType
  f::F
  uprev::uType
  kprev::ksEltype
  tprev::tType
  Ts::Vector{tType}
  tableau::ExplicitRKTableau
  autodiff::Bool
  adaptiveorder::Int
  order::Int
  fsal::Bool
  alg::algType
  custom_callback::Bool
  rate_prototype::rateType
  notsaveat_idxs::Vector{Int}
  calcprevs::Bool
  dtcache::tType
  dt_mod::tType
  iter::Int
  saveiter::Int
  saveiter_dense::Int
  cursaveat::Int
  event_cache::ECType
  kshortsize::Int
  reeval_fsal::Bool
  opts::O
end

@def ode_preamble begin
  local t::tType
  local dt::tType

  @unpack u,k,t,dt,Ts,autodiff,fsal,alg,rate_prototype = integrator
  f = integrator.f

  sizeu = size(u)
  Tfinal = Ts[end]

  local T::tType
  local utmp::uType

  # Setup FSAL
  if uType <: Number
    utmp = zero(uType)
    fsallast = zero(rateType)
    if fsal
      fsalfirst = f(t,u)
    else
      fsalfirst = zero(rateType)
    end
  else
    utmp = zeros(u)
    fsallast = similar(rate_prototype)
    fsalfirst = similar(rate_prototype)
    if fsal
      f(t,u,fsalfirst)
    end
  end

  tTypeNoUnits = typeof(integrator.opts.qoldinit)
  uEltypeNoUnits = typeof(integrator.opts.reltol)
  uEltype = typeof(integrator.opts.abstol)
  local Θ = one(t)/one(t) # No units
  local q::tTypeNoUnits = 0
  local dtpropose::tType = tType(0)
  local q11::tTypeNoUnits = 0
  local qold::tTypeNoUnits = integrator.opts.qoldinit
  if uType <: Number || !integrator.custom_callback
    cache = ()
  end
  qminc = inv(integrator.opts.qmin) #facc1
  qmaxc = inv(integrator.opts.qmax) #facc2
  local EEst::tTypeNoUnits = zero(t)
  integrator.opts.progress && (prog = Juno.ProgressBar(name=integrator.opts.progress_name))
end

@def ode_loopheader begin
  integrator.iter += 1

  if integrator.opts.adaptive
    dt = min(dt,abs(T-t)) # Step to the end
  elseif integrator.dtcache == 0 # Use tstops
    dt = abs(T-t)
  else # always try to step with dtcache
    dt = min(integrator.dtcache,abs(T-t)) # Step to the end
  end

  if integrator.iter > integrator.opts.maxiters
    warn("Interrupted. Larger maxiters is needed.")
    @ode_postamble
  end
  if dt == tType(0)
    warn("dt == 0. Aborting")
    @ode_postamble
  end
  if any(isnan,u)
    warn("NaNs detected. Aborting")
    @ode_postamble
  end

  if uType<:AbstractArray && integrator.custom_callback
    uidx = eachindex(u)
  end
end

@def ode_savevalues begin
  if !isempty(integrator.opts.saveat) # Perform saveat
    while integrator.cursaveat <= length(integrator.opts.saveat) && integrator.opts.saveat[integrator.cursaveat]<= t
      integrator.saveiter += 1
      if integrator.opts.saveat[integrator.cursaveat]<t # If <t, interpolate
        curt = integrator.opts.saveat[integrator.cursaveat]
        ode_addsteps!(integrator.k,integrator.tprev,integrator.uprev,integrator.dt,integrator.alg,integrator.f)
        Θ = (curt - integrator.tprev)/integrator.dt
        val = ode_interpolant(Θ,integrator.dt,integrator.uprev,u,integrator.kprev,integrator.k,integrator.alg)
        copyat_or_push!(integrator.sol.t,integrator.saveiter,curt)
        copyat_or_push!(integrator.sol.u,integrator.saveiter,val)
      else # ==t, just save
        copyat_or_push!(integrator.sol.t,integrator.saveiter,t)
        copyat_or_push!(integrator.sol.u,integrator.saveiter,u)
        if integrator.opts.dense
          integrator.saveiter_dense += 1
          copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
          copyat_or_push!(integrator.notsaveat_idxs,integrator.saveiter_dense,integrator.saveiter)
        end
      end
      integrator.cursaveat+=1
    end
  end
  if integrator.opts.save_timeseries && integrator.iter%integrator.opts.timeseries_steps==0
    integrator.saveiter += 1
    copyat_or_push!(integrator.sol.u,integrator.saveiter,u)
    copyat_or_push!(integrator.sol.t,integrator.saveiter,t)
    if integrator.opts.dense
      integrator.saveiter_dense += 1
      copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
      copyat_or_push!(integrator.notsaveat_idxs,integrator.saveiter_dense,integrator.saveiter)
    end
  end
  if isspecialdense(alg)
    resize!(integrator.k,integrator.kshortsize)
  end
end

@def ode_postamble begin
  if integrator.sol.t[end] != t
    integrator.saveiter += 1
    copyat_or_push!(integrator.sol.t,integrator.saveiter,t)
    copyat_or_push!(integrator.sol.u,integrator.saveiter,u)
    if integrator.opts.dense
      integrator.saveiter_dense +=1
      copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
      copyat_or_push!(integrator.notsaveat_idxs,integrator.saveiter_dense,integrator.saveiter)
    end
  end
  integrator.opts.progress && Juno.done(prog)
end

@def pack_integrator begin
  integrator.dt = dt
  integrator.dt_mod = tType(1)
  integrator.k = k
end

@def unpack_integrator begin
end

@def ode_loopfooter begin
  if integrator.opts.adaptive
    q11 = EEst^integrator.opts.beta1
    q = q11/(qold^integrator.opts.beta2)
    q = max(qmaxc,min(qminc,q/integrator.opts.gamma))
    dtnew = dt/q
    ttmp = t + dt
    if !integrator.opts.isoutofdomain(ttmp,utmp) && EEst <= 1.0 # Accept
      t = ttmp
      if uType <: AbstractArray # Treat mutables differently
        recursivecopy!(u, utmp)
      else
        u = utmp
      end

      qold = max(EEst,integrator.opts.qoldinit)
      dtpropose = min(integrator.opts.dtmax,dtnew)
      @pack_integrator
      if integrator.custom_callback
        t,T = integrator.opts.callback(alg,f,t,u,dt,cache,T,Ts,integrator)
      else
        @ode_savevalues
      end
      @unpack_integrator
      dt = integrator.dt_mod*max(dtpropose,integrator.opts.dtmin) #abs to fix complex sqrt issue at end

      if fsal
        if integrator.reeval_fsal
          if uType <: AbstractArray
            f(t,u,fsalfirst)
          else
            fsalfirst = f(t,u)
          end
          integrator.reeval_fsal = false
        else
          if uType <: AbstractArray
            recursivecopy!(fsalfirst,fsallast)
          else
            fsalfirst = fsallast
          end
        end
      end

      if integrator.calcprevs
        # Store previous for interpolation
        integrator.tprev = t
        if uType <: AbstractArray
          recursivecopy!(integrator.uprev,u)
        else
          integrator.uprev = u
        end
        if integrator.opts.calck
          if ksEltype <: AbstractArray
            recursivecopy!(integrator.kprev,k)
          else
            integrator.kprev = k
          end
        end
      end
    else # Reject
      dt = dt/min(qminc,q11/integrator.opts.gamma)
    end
  else #Not adaptive
    t += dt
    if fsal
      if uType <: AbstractArray
        recursivecopy!(fsalfirst,fsallast)
      else
        fsalfirst = fsallast
      end
    end
    @pack_integrator
    if integrator.custom_callback
      t,T = integrator.opts.callback(alg,f,t,u,dt,cache,T,Ts,integrator)
    else
      @ode_savevalues
    end
    @unpack_integrator
    dt *= integrator.dt_mod
    if integrator.calcprevs
      # Store previous for interpolation
      integrator.tprev = t
      if uType <: AbstractArray
        recursivecopy!(integrator.uprev,u)
      else
        integrator.uprev = u
      end
      if integrator.opts.calck
        if ksEltype <: AbstractArray
          recursivecopy!(integrator.kprev,k)
        else
          integrator.kprev = k
        end
      end
    end
  end
  if integrator.opts.progress && integrator.iter%integrator.opts.progress_steps==0
    Juno.msg(prog,integrator.opts.progress_message(dt,t,u))
    Juno.progress(prog,t/Tfinal)
  end
end
