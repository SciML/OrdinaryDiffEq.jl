function solve{uType,tType,isinplace,T<:OrdinaryDiffEqAlgorithm,F}(
  prob::AbstractODEProblem{uType,tType,Val{isinplace},F},
  algType::Type{T},timeseries=[],ts=[],ks=[];
  dt = 0.0,save_timeseries = true,
  timeseries_steps = 1,tableau = ODE_DEFAULT_TABLEAU,
  dense = true,calck = nothing,alg_hint = :nonstiff,
  timeseries_errors = true,dense_errors = false,
  saveat = tType[],tstops = tType[],
  adaptive = true,gamma=.9,abstol=1//10^6,reltol=1//10^3,
  qmax=nothing,qmin=nothing,qoldinit=1//10^4, fullnormalize=true,
  beta2=nothing,beta1=nothing,maxiters = 10000,
  dtmax=tType((prob.tspan[end]-prob.tspan[1])),
  dtmin=tType <: AbstractFloat ? tType(10)*eps(tType) : tType(1//10^(10)),
  autodiff=true,internalnorm = ODE_DEFAULT_NORM,
  progressbar=false,progress_steps=1000,progressbar_name="ODE",
  callback=nothing,kwargs...)

  alg = algType()

  tspan = prob.tspan

  if tspan[end]-tspan[1]<tType(0)
    error("final time must be greater than starting time. Aborting.")
  end

  t = tspan[1]
  Ts = sort(unique([tstops;tspan[2]]))

  if tspan[end] < Ts[end]
      error("Final saving timepoint is past the solving timespan")
  end
  if t > Ts[1]
      error("First saving timepoint is before the solving timespan")
  end


  atomloaded = isdefined(Main,:Atom)
  u0 = prob.u0
  uEltype = eltype(u0)

  # Get the control variables

  if callback == nothing
    callback = ODE_DEFAULT_CALLBACK
    custom_callback = false
  else
    custom_callback = true
  end

  if uEltype<:Number
    u = copy(u0)
  else
    u = deepcopy(u0)
  end

  ks = Vector{uType}(0)

  order = alg.order
  adaptiveorder = 0

  if typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm
    adaptiveorder = alg.adaptiveorder
    if adaptive == true
      dt = 1.0*dt # Convert to float in a way that keeps units
    end
  else
    adaptive = false
  end

  if typeof(alg) <: ExplicitRK
    @unpack order,adaptiveorder = tableau
  end

  if !isinplace && typeof(u)<:AbstractArray
    f! = (t,u,du) -> (du[:] = prob.f(t,u))
  else
    f! = prob.f
  end

  uEltypeNoUnits = typeof(recursive_one(u))
  tTypeNoUnits   = typeof(recursive_one(t))

  if dt==0
    dt = ode_determine_initdt(u0,t,uEltype(abstol),uEltypeNoUnits(reltol),internalnorm,f!,order)
  end

  rate_prototype = u/zero(t)
  rateType = typeof(rate_prototype) ## Can be different if united

  saveat = tType[convert(tType,x) for x in setdiff(saveat,tspan)]

  if calck==nothing
    calck = !isempty(saveat) || dense
  end

  ### Algorithm-specific defaults ###

  if qmin == nothing # Use default qmin
    if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
      qmin = 0.2
    elseif typeof(alg) <: DP8
      qmin = 0.333
    else
      qmin = 0.2
    end
  end
  if qmax == nothing # Use default qmax
    if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
      qmax = 10.0
    elseif typeof(alg) <: DP8
      qmax = 6.0
    else
      qmax = 10.0
    end
  end
  if beta2 == nothing # Use default β₂
    if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
      β₂ = 0.04
    elseif typeof(alg) <: DP8
      β₂ = 0.00
    else
      β₂ = 0.4 / order
    end
  else
    β₂ = beta2
  end
  if beta1 == nothing # Use default β₁
    if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
      β₁ = 1/order - .75β₂
    elseif typeof(alg) <: DP8
      β₁ = 1/order - .2β₂
    else
      β₁ = .7/order
    end
  else
    β₁ = beta1
  end

  fsal = false
  if isfsal(alg)
    fsal = true
  elseif typeof(alg) <: ExplicitRK
    @unpack fsal = tableau
  end

  abstol = uEltype(1)*abstol

  if isspecialdense(alg)
    ksEltype = Vector{rateType} # Store more ks for the special algs
  else
    ksEltype = rateType # Makes simple_dense
  end

  timeseries = convert(Vector{uType},timeseries)
  ts = convert(Vector{tType},ts)
  ks = convert(Vector{ksEltype},ks)
  if length(timeseries) == 0
    push!(timeseries,copy(u))
  else
    timeseries[1] = copy(u)
  end

  if length(ts) == 0
    push!(ts,t)
  else
    timeseries[1] = copy(u)
  end

  if ksEltype == rateType
    if uType <: Number
      rate_prototype = f!(t,u)
    else
      f!(t,u,rate_prototype)
    end
    push!(ks,rate_prototype)
  else # Just push a dummy in for special dense since first is not used.
    push!(ks,[rate_prototype])
  end
  γ = gamma
  # @code_warntype ode_solve(ODEIntegrator{alg,uType,uEltype,ndims(u)+1,tType,tTypeNoUnits,uEltypeNoUnits,rateType,ksEltype}(timeseries,ts,ks,f!,u,t,k,dt,Ts,maxiters,timeseries_steps,save_timeseries,adaptive,abstol,reltol,γ,qmax,qmin,dtmax,dtmin,internalnorm,progressbar,tableau,autodiff,adaptiveorder,order,atomloaded,progress_steps,β₁,β₂,qoldinit,fsal,dense,saveat,alg,callback,custom_callback,calck))
  u,t = ode_solve(ODEIntegrator{typeof(alg),uType,uEltype,ndims(u)+1,tType,uEltypeNoUnits,tTypeNoUnits,rateType,ksEltype,typeof(f!),typeof(internalnorm),typeof(callback)}(timeseries,ts,ks,f!,u,t,dt,Ts,maxiters,timeseries_steps,save_timeseries,adaptive,abstol,reltol,γ,qmax,qmin,dtmax,dtmin,internalnorm,progressbar,tableau,autodiff,adaptiveorder,order,atomloaded,progress_steps,progressbar_name,β₁,β₂,qoldinit,fsal,dense,saveat,alg,callback,custom_callback,calck))


  saveat_idxs = find((x)->x∈saveat,ts)
  t_nosaveat = view(ts,symdiff(1:length(ts),saveat_idxs))
  u_nosaveat = view(timeseries,symdiff(1:length(ts),saveat_idxs))

  if dense
    interp = (tvals) -> ode_interpolation(tvals,t_nosaveat,u_nosaveat,ks,alg,f!)
  else
    interp = (tvals) -> nothing
  end

  build_ode_solution(prob,alg,ts,timeseries,
                    dense=dense,k=ks,interp=interp,
                    timeseries_errors = timeseries_errors,
                    dense_errors = dense_errors)
end

function ode_determine_initdt{uType,tType,uEltypeNoUnits}(u0::uType,t::tType,abstol,reltol::uEltypeNoUnits,internalnorm,f,order)
  f₀ = similar(u0./t); f₁ = similar(u0./t); u₁ = similar(u0)
  d₀ = internalnorm(u0./(abstol+u0*reltol))
  f(t,u0,f₀)
  d₁ = internalnorm(f₀./(abstol+u0*reltol)*tType(1))/tType(1)
  T0 = typeof(d₀)
  T1 = typeof(d₁)
  if d₀ < T0(1//10^(5)) || d₁ < T1(1//10^(5))
    dt₀ = tType(1//10^(6))
  else
    dt₀ = tType((d₀/d₁)/100)
  end
  @inbounds for i in eachindex(u0)
     u₁[i] = u0[i] + dt₀*f₀[i]
  end
  f(t+dt₀,u₁,f₁)
  d₂ = internalnorm((f₁.-f₀)./(abstol+u0*reltol)*tType(1))/dt₀
  if max(d₁,d₂)<=T1(1//10^(15))
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(max(d₁,d₂)/T1(1)))/(order)))
  end
  dt = min(100dt₀,dt₁)
end

function ode_determine_initdt{uType<:Number,tType,uEltypeNoUnits}(u0::uType,t::tType,abstol,reltol::uEltypeNoUnits,internalnorm,f,order)
  d₀ = abs(u0./(abstol+u0*reltol))
  f₀ = f(t,u0)
  d₁ = abs(f₀./(abstol+u0*reltol))
  T0 = typeof(d₀)
  T1 = typeof(d₁)
  if d₀ < T0(1//10^(5)) || d₁ < T1(1//10^(5))
    dt₀ = tType(1//10^(6))
  else
    dt₀ = tType((d₀/d₁)/100)
  end
  u₁ = u0 + dt₀*f₀
  f₁ = f(t+dt₀,u₁)
  d₂ = abs((f₁-f₀)./(abstol+u0*reltol))/dt₀*tType(1)
  if max(d₁,d₂) <= T1(1//10^(15))
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(max(d₁,d₂)/T1(1)))/(order)))
  end
  dt = min(100dt₀,dt₁)
end
