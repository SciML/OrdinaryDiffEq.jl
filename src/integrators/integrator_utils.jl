immutable ODEIntegrator{algType<:OrdinaryDiffEqAlgorithm,uType<:Union{AbstractArray,Number},uEltype,N,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,F,F2,F3,F4}
  timeseries::Vector{uType}
  ts::Vector{tType}
  ks::Vector{ksEltype}
  f::F
  u::uType
  t::tType
  dt::tType
  Ts::Vector{tType}
  maxiters::Int
  timeseries_steps::Int
  save_timeseries::Bool
  adaptive::Bool
  abstol::uEltype
  reltol::uEltypeNoUnits
  γ::tTypeNoUnits
  qmax::tTypeNoUnits
  qmin::tTypeNoUnits
  dtmax::tType
  dtmin::tType
  internalnorm::F2
  progressbar::Bool
  tableau::ExplicitRKTableau
  autodiff::Bool
  adaptiveorder::Int
  order::Int
  atomloaded::Bool
  progress_steps::Int
  progressbar_name::String
  β₁::tTypeNoUnits
  β₂::tTypeNoUnits
  qoldinit::tTypeNoUnits
  fsal::Bool
  dense::Bool
  saveat::Vector{tType}
  alg::algType
  callback::F3
  isoutofdomain::F4
  custom_callback::Bool
  calck::Bool
end

@def ode_preamble begin
  local u::uType
  local t::tType
  local dt::tType
  local Ts::Vector{tType}
  local adaptiveorder::Int
  @unpack f,u,t,dt,Ts,maxiters,timeseries_steps,γ,qmax,qmin,save_timeseries,adaptive,progressbar,autodiff,adaptiveorder,order,atomloaded,progress_steps,progressbar_name,β₂,β₁,qoldinit,fsal, dense, saveat, alg, callback, isoutofdomain, custom_callback,calck = integrator
  timeseries = integrator.timeseries
  ts = integrator.ts
  ks = integrator.ks
  const calcprevs = !isempty(saveat) || custom_callback # Calculate the previous values
  const issimple_dense = (ksEltype==rateType) # Means ks[i] = f(t[i],timeseries[i]), for Hermite

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
  local rate_prototype::rateType
  if ksEltype <: AbstractArray  && !issimple_dense
    k = ksEltype[]
    kprev = ksEltype[]
    rate_prototype = ks[1][1]
  elseif ksEltype <: Number
    k = ksEltype(0)
    kprev = ksEltype(0)
    rate_prototype = ks[1][1]
  else # it is simple_dense
    k = ksEltype(zeros(Int64,N-1)...) # Needs the zero for dimension 3+
    kprev = ksEltype(zeros(Int64,N-1)...)
    rate_prototype = ks[1]
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

  local cursaveat::Int = 1
  local Θ = one(t)/one(t) # No units
  local tprev::tType = t
  local uprev::uType = deepcopy(u)
  local standard::uEltype = uEltype(0)
  local q::tTypeNoUnits = 0
  local dtpropose::tType = tType(0)
  local q11::tTypeNoUnits = 0
  local qold::tTypeNoUnits = qoldinit
  local β₁::tTypeNoUnits
  if uType <: Number || !custom_callback
    cache = ()
  end
  qminc = inv(qmin) #facc1
  qmaxc = inv(qmax) #facc2
  local EEst::tTypeNoUnits = zero(uEltypeNoUnits)
  if adaptive
    @unpack abstol,reltol,qmax,dtmax,dtmin,internalnorm = integrator
  end

  if issimple_dense #If issimple_dense, then ks[1]=f(ts[1],timeseries[1])
    const kshortsize = 1
    if calck
      if ksEltype <: AbstractArray
        k = similar(rate_prototype)
      end
      kprev = copy(k)
    end
  end ## if not simple_dense, you have to initialize k and push the ks[1]!

  progressbar && (prog = ProgressBar(name=progressbar_name))
end

@def ode_loopheader begin
  iter += 1
  if iter > maxiters
    warn("Interrupted. Larger maxiters is needed.")
    @ode_postamble
  end
  dt = min(dt,abs(T-t))
  if uType<:AbstractArray && custom_callback
    uidx = eachindex(u)
  end
end

@def ode_savevalues begin
  if !isempty(saveat) # Perform saveat
    while cursaveat <= length(saveat) && saveat[cursaveat]<= t
      if saveat[cursaveat]<t # If we already saved at the point, ignore it
        saveiter += 1
        curt = saveat[cursaveat]
        ode_addsteps!(k,tprev,uprev,dt,alg,f)
        Θ = (curt - tprev)/dt
        val = ode_interpolant(Θ,dt,uprev,u,kprev,k,alg)
        copyat_or_push!(ts,saveiter,curt)
        copyat_or_push!(timeseries,saveiter,val)
      end
      cursaveat+=1
    end
  end
  if save_timeseries && iter%timeseries_steps==0
    saveiter += 1
    copyat_or_push!(timeseries,saveiter,u)
    copyat_or_push!(ts,saveiter,t)
    if dense
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
  progressbar && done(prog)
  return u,t
end

@def ode_loopfooter begin
  if adaptive
    q11 = EEst^β₁
    q = q11/(qold^β₂)
    q = max(qmaxc,min(qminc,q/γ))
    dtnew = dt/q
    ttmp = t + dt
    if !isoutofdomain(ttmp,utmp) && EEst <= 1.0 # Accept
      t = ttmp
      if uType <: AbstractArray # Treat mutables differently
        recursivecopy!(u, utmp)
      else
        u = utmp
      end

      qold = max(EEst,qoldinit)
      dtpropose = min(dtmax,dtnew)
      dtprev = dt
      dt = max(dtpropose,dtmin) #abs to fix complex sqrt issue at end
      cursaveat,saveiter,dt,t,T,reeval_fsal = callback(alg,f,t,u,k,tprev,uprev,kprev,ts,timeseries,ks,dtprev,dt,saveat,cursaveat,saveiter,iter,save_timeseries,timeseries_steps,uEltype,ksEltype,dense,kshortsize,issimple_dense,fsal,fsalfirst,cache,calck,T,Ts)

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
        if calck
          if ksEltype <: AbstractArray
            recursivecopy!(kprev,k)
          else
            kprev = k
          end
        end
      end
    else # Reject
      dt = dt/min(qminc,q11/γ)
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
    cursaveat,saveiter,dt,t,T,reeval_fsal = callback(alg,f,t,u,k,tprev,uprev,kprev,ts,timeseries,ks,dtprev,dt,saveat,cursaveat,saveiter,iter,save_timeseries,timeseries_steps,uEltype,ksEltype,dense,kshortsize,issimple_dense,fsal,fsalfirst,cache,calck,T,Ts)
    if calcprevs
      # Store previous for interpolation
      tprev = t
      if uType <: AbstractArray
        recursivecopy!(uprev,u)
      else
        uprev = u
      end
      if calck
        if ksEltype <: AbstractArray
          recursivecopy!(kprev,k)
        else
          kprev = k
        end
      end
    end
  end
  if progressbar && iter%progress_steps==0
    msg(prog,"dt="*string(dt))
    progress(prog,t/Tfinal)
  end
end
