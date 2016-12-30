function solve{uType,tType,isinplace,algType<:OrdinaryDiffEqAlgorithm,F}(
  prob::AbstractODEProblem{uType,tType,isinplace,F},
  alg::algType,timeseries=[],ts=[],ks=[];
  timeseries_errors = true,dense_errors = false,
  kwargs...)

  integrator = init(prob,alg,timeseries,ts,ks;kwargs...)
  solve!(integrator,timeseries_errors=timeseries_errors,dense_errors=dense_errors)
  integrator.sol
end

function init{uType,tType,isinplace,algType<:OrdinaryDiffEqAlgorithm,F}(
  prob::AbstractODEProblem{uType,tType,isinplace,F},
  alg::algType,timeseries=[],ts=[],ks=[];
  dt = tType(0),save_timeseries = true,
  timeseries_steps = 1,tableau = ODE_DEFAULT_TABLEAU,
  dense = save_timeseries,calck = nothing,alg_hint = :nonstiff,
  saveat = tType[],tstops = tType[],
  adaptive = true,
  gamma=.9,
  abstol=1//10^6,
  reltol=1//10^3,
  qmax=qmax_default(alg),qmin=qmin_default(alg),
  qoldinit=1//10^4, fullnormalize=true,
  beta2=beta2_default(alg),
  beta1=beta1_default(alg,beta2),
  maxiters = 1000000,
  dtmax=tType((prob.tspan[end]-prob.tspan[1])),
  dtmin=tType <: AbstractFloat ? tType(10)*eps(tType) : tType(1//10^(10)),
  autodiff=false,internalnorm = ODE_DEFAULT_NORM,
  isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
  progress=false,progress_steps=1000,progress_name="ODE",
  progress_message = ODE_DEFAULT_PROG_MESSAGE,
  event_cache=nothing,callback=nothing,kwargs...)

  tspan = prob.tspan

  if tspan[end]-tspan[1]<tType(0)
    error("Final time must be greater than starting time. Aborting.")
  end

  t = tspan[1]
  Ts = sort(unique([tstops;tspan[2]]))

  if tspan[end] < Ts[end]
      error("Final saving timepoint is past the solving timespan")
  end
  if t > Ts[1]
      error("First saving timepoint is before the solving timespan")
  end

  if !(typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm) && dt == 0 && isempty(tstops)
      error("Fixed timestep methods require a choice of dt or choosing the tstops")
  end

  u0 = prob.u0
  uEltype = eltype(u0)

  # Get the control variables

  if callback == nothing
    callback = ODE_DEFAULT_CALLBACK
    custom_callback = false
  else
    custom_callback = true
  end

  if uType<:Array || uType <: Number
    u = copy(u0)
  else
    u = deepcopy(u0)
  end

  ks = Vector{uType}(0)

  order = alg.order
  adaptiveorder = 0

  if typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm
    adaptiveorder = alg.adaptiveorder
    #=
    if adaptive == true
      dt = 1.0*dt # Convert to float in a way that keeps units
    end
    =#
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

  if dt==0 && adaptive
    dt = tType(ode_determine_initdt(u0,t,uEltype(abstol),uEltypeNoUnits(reltol),internalnorm,f!,order))
  end

  rate_prototype = u/zero(t)
  rateType = typeof(rate_prototype) ## Can be different if united

  saveat = tType[convert(tType,x) for x in setdiff(saveat,tspan)]

  if calck==nothing
    calck = !isempty(saveat) || dense
  end

  ### Algorithm-specific defaults ###

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

  opts = DEOptions(maxiters,timeseries_steps,save_timeseries,adaptive,uEltype(abstol),
    uEltypeNoUnits(reltol),gamma,qmax,qmin,dtmax,dtmin,internalnorm,progress,progress_steps,
    progress_name,progress_message,beta1,beta2,tTypeNoUnits(qoldinit),dense,saveat,
    callback,isoutofdomain,calck)

  notsaveat_idxs = Int[1]

  if ksEltype <: AbstractArray  &&  isspecialdense(alg)
    k = ksEltype[]
    kprev = ksEltype[]
  elseif ksEltype <: Number
    k = ksEltype(0)
    kprev = ksEltype(0)
  else # it is simple_dense
    k = ksEltype(zeros(Int64,ndims(u))...) # Needs the zero for dimension 3+
    kprev = ksEltype(zeros(Int64,ndims(u))...)
  end

  if !isspecialdense(alg) #If issimple_dense, then ks[1]=f(ts[1],timeseries[1])
    if calck
      if ksEltype <: AbstractArray
        k = similar(rate_prototype)
      end
      kprev = copy(k)
    end
  end ## if not simple_dense, you have to initialize k and push the ks[1]!

  if uType <: Array
    uprev = copy(u)
  else
    uprev = deepcopy(u)
  end

  if dense
    #notsaveat_idxs  = find((x)->(x∉saveat)||(x∈Ts),ts)
    id = InterpolationData(f!,timeseries,ts,ks,notsaveat_idxs)
    interp = (tvals) -> ode_interpolation(alg,tvals,id)
  else
    interp = (tvals) -> nothing
  end

  sol = build_solution(prob,alg,ts,timeseries,
                    dense=dense,k=ks,interp=interp,
                    calculate_error = false)

  calcprevs = !isempty(saveat) || custom_callback # Calculate the previous values
  tprev = t
  dtcache = tType(dt)
  dt_mod = tType(1)
  iter = 0
  saveiter = 1 # Starts at 1 so first save is at 2
  saveiter_dense = 1
  cursaveat = 1
  kshortsize = 1
  reeval_fsal = false

  integrator = ODEIntegrator{algType,uType,tType,eltype(ks),typeof(tableau),typeof(sol),
                             typeof(rate_prototype),typeof(f!),typeof(event_cache),typeof(opts)}(
                             sol,u,k,t,tType(dt),f!,uprev,kprev,tprev,
                             Ts,tableau,autodiff,adaptiveorder,order,fsal,
                             alg,custom_callback,rate_prototype,
                             notsaveat_idxs,calcprevs,dtcache,dt_mod,
                             iter,saveiter,saveiter_dense,cursaveat,
                             event_cache,kshortsize,reeval_fsal,opts)
  integrator
end

function solve!(integrator::ODEIntegrator;timeseries_errors = true,dense_errors = false)
  ode_solve(integrator)

  if typeof(integrator.sol.prob) <: AbstractODETestProblem
    calculate_solution_errors!(integrator.sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
  end
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
