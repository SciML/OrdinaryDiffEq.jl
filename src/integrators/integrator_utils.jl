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

immutable ODEIntegrator{algType<:OrdinaryDiffEqAlgorithm,uType<:Union{AbstractArray,Number},tType,ksEltype,F,rateType,O}
  timeseries::Vector{uType}
  ts::Vector{tType}
  ks::Vector{ksEltype}
  f::F
  u::uType
  t::tType
  dt::tType
  Ts::Vector{tType}
  tableau::ExplicitRKTableau
  autodiff::Bool
  adaptiveorder::Int
  order::Int
  fsal::Bool
  alg::algType
  custom_callback::Bool
  rate_prototype::rateType
  opts::O
end

@def ode_preamble begin
  local u::uType
  local t::tType
  local dt::tType
  local Ts::Vector{tType}
  local adaptiveorder::Int
  @unpack f,u,t,dt,Ts,autodiff,adaptiveorder,order,fsal,alg,custom_callback,rate_prototype = integrator
  timeseries = integrator.timeseries
  ts = integrator.ts
  ks = integrator.ks
  const calcprevs = !isempty(integrator.opts.saveat) || custom_callback # Calculate the previous values
  const issimple_dense = !isspecialdense(alg)
  const dtcache = dt
  # Need to initiate ks in the method



  Tfinal = Ts[end]
  local iter::Int = 0
  local saveiter::Int = 1 # Starts at 1 so first save is at 2
  local T::tType
  sizeu = size(u)
  local utmp::uType
  local k::ksEltype
  local kprev::ksEltype
  local kshortsize::Int
  local reeval_fsal::Bool
  if ksEltype <: AbstractArray  && !issimple_dense
    k = ksEltype[]
    kprev = ksEltype[]
  elseif ksEltype <: Number
    k = ksEltype(0)
    kprev = ksEltype(0)
  else # it is simple_dense
    k = ksEltype(zeros(Int64,ndims(u))...) # Needs the zero for dimension 3+
    kprev = ksEltype(zeros(Int64,ndims(u))...)
  end

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
  local cursaveat::Int = 1
  local Θ = one(t)/one(t) # No units
  local tprev::tType = t
  local uprev::uType = deepcopy(u)
  local q::tTypeNoUnits = 0
  local dtpropose::tType = tType(0)
  local q11::tTypeNoUnits = 0
  local qold::tTypeNoUnits = integrator.opts.qoldinit
  if uType <: Number || !custom_callback
    cache = ()
  end
  qminc = inv(integrator.opts.qmin) #facc1
  qmaxc = inv(integrator.opts.qmax) #facc2
  local EEst::tTypeNoUnits = zero(uEltypeNoUnits)

  if issimple_dense #If issimple_dense, then ks[1]=f(ts[1],timeseries[1])
    const kshortsize = 1
    if integrator.opts.calck
      if ksEltype <: AbstractArray
        k = similar(rate_prototype)
      end
      kprev = copy(k)
    end
  end ## if not simple_dense, you have to initialize k and push the ks[1]!

  integrator.opts.progress && (prog = Juno.ProgressBar(name=integrator.opts.progress_name))
end

@def ode_loopheader begin
  iter += 1

  if integrator.opts.adaptive
    dt = min(dt,abs(T-t)) # Step to the end
  elseif dtcache == 0 # Use tstops
    dt = abs(T-t)
  else # always try to step with dtcache
    dt = min(dtcache,abs(T-t)) # Step to the end
  end

  if iter > integrator.opts.maxiters
    warn("Interrupted. Larger maxiters is needed.")
    @ode_postamble
  end
  if dt == 0
    warn("dt == 0. Aborting")
    @ode_postamble
  end
  if any(isnan,u)
    warn("NaNs detected. Aborting")
    @ode_postamble
  end

  if uType<:AbstractArray && custom_callback
    uidx = eachindex(u)
  end
end

@def ode_savevalues begin
  if !isempty(integrator.opts.saveat) # Perform saveat
    while cursaveat <= length(integrator.opts.saveat) && integrator.opts.saveat[cursaveat]<= t
      saveiter += 1
      if integrator.opts.saveat[cursaveat]<t # If <t, interpolate
        curt = integrator.opts.saveat[cursaveat]
        ode_addsteps!(k,tprev,uprev,dtprev,alg,f)
        Θ = (curt - tprev)/dtprev
        val = ode_interpolant(Θ,dtprev,uprev,u,kprev,k,alg)
        copyat_or_push!(ts,saveiter,curt)
        copyat_or_push!(timeseries,saveiter,val)
      else # ==t, just save
        copyat_or_push!(ts,saveiter,t)
        copyat_or_push!(timeseries,saveiter,u)
        if integrator.opts.dense
          copyat_or_push!(ks,saveiter,k)
        end
      end
      cursaveat+=1
    end
  end
  if integrator.opts.save_timeseries && iter%integrator.opts.timeseries_steps==0
    saveiter += 1
    copyat_or_push!(timeseries,saveiter,u)
    copyat_or_push!(ts,saveiter,t)
    if integrator.opts.dense
      copyat_or_push!(ks,saveiter,k)
    end
  end
  if !issimple_dense
    resize!(k,kshortsize)
  end
end

@def ode_postamble begin
  if ts[end] != t
    push!(ts,t)
    push!(timeseries,u)
  end
  integrator.opts.progress && Juno.done(prog)
  return u,t
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
      dtprev = dt
      dt = max(dtpropose,integrator.opts.dtmin) #abs to fix complex sqrt issue at end
      if custom_callback
        cursaveat,saveiter,dt,t,T,reeval_fsal = integrator.opts.callback(alg,f,t,u,k,tprev,uprev,kprev,ts,timeseries,ks,dtprev,dt,cursaveat,saveiter,iter,uEltype,ksEltype,kshortsize,issimple_dense,fsal,fsalfirst,cache,T,Ts,integrator)
      else
        @ode_savevalues
        reeval_fsal = false
      end

      if fsal
        if reeval_fsal
          if uType <: AbstractArray
            f(t,u,fsalfirst)
          else
            fsalfirst = f(t,u)
          end
        else
          if uType <: AbstractArray
            recursivecopy!(fsalfirst,fsallast)
          else
            fsalfirst = fsallast
          end
        end
      end

      if calcprevs
        # Store previous for interpolation
        tprev = t
        if uType <: AbstractArray
          recursivecopy!(uprev,u)
        else
          uprev = u
        end
        if integrator.opts.calck
          if ksEltype <: AbstractArray
            recursivecopy!(kprev,k)
          else
            kprev = k
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
    dtprev = dt
    if custom_callback
      cursaveat,saveiter,dt,t,T,reeval_fsal = integrator.opts.callback(alg,f,t,u,k,tprev,uprev,kprev,ts,timeseries,ks,dtprev,dt,cursaveat,saveiter,iter,uEltype,ksEltype,kshortsize,issimple_dense,fsal,fsalfirst,cache,T,Ts,integrator)
    else
      @ode_savevalues
      reeval_fsal = false
    end
    if calcprevs
      # Store previous for interpolation
      tprev = t
      if uType <: AbstractArray
        recursivecopy!(uprev,u)
      else
        uprev = u
      end
      if integrator.opts.calck
        if ksEltype <: AbstractArray
          recursivecopy!(kprev,k)
        else
          kprev = k
        end
      end
    end
  end
  if integrator.opts.progress && iter%integrator.opts.progress_steps==0
    Juno.msg(prog,integrator.opts.progress_message(dt,t,u))
    Juno.progress(prog,t/Tfinal)
  end
end
