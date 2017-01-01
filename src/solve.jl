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
  alg::algType,timeseries_init=uType[],ts_init=tType[],ks_init=[];
  dt = tType(0),save_timeseries = true,
  timeseries_steps = 1,
  dense = save_timeseries,
  saveat = tType[],tstops = tType[],
  calck = (!isempty(setdiff(saveat,tstops)) || dense),
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
  internalnorm = ODE_DEFAULT_NORM,
  isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
  progress=false,progress_steps=1000,progress_name="ODE",
  progress_message = ODE_DEFAULT_PROG_MESSAGE,
  event_cache=nothing,callback=nothing,kwargs...)

  tspan = prob.tspan
  tdir = sign(tspan[end]-tspan[1])

  t = tspan[1]

  if !(typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm) && dt == tType(0) && isempty(tstops)
      error("Fixed timestep methods require a choice of dt or choosing the tstops")
  end

  Ts = push!(tstops,tspan[2])

  #Heapify first
  if tdir*tspan[end] < tdir*Ts[end]
      error("Final saving timepoint is past the solving timespan")
  end
  if tdir*t > tdir*Ts[1]
      error("First saving timepoint is before the solving timespan")
  end

  u0 = prob.u0
  uEltype = eltype(u0)

  # Get the control variables

  (uType<:Array || uType <: Number) ? u = copy(u0) : u = deepcopy(u0)

  ks = Vector{uType}(0)

  order = alg_order(alg)
  adaptiveorder = 0

  typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm ? adaptiveorder = alg_adaptive_order(alg) : adaptive = false

  if typeof(alg) <: ExplicitRK
    @unpack order,adaptiveorder = alg.tableau
  end

  if !isinplace && typeof(u)<:AbstractArray
    f! = (t,u,du) -> (du[:] = prob.f(t,u))
  else
    f! = prob.f
  end

  uEltypeNoUnits = typeof(recursive_one(u))
  tTypeNoUnits   = typeof(recursive_one(t))

  if dt == zero(dt) && adaptive
    dt = tType(ode_determine_initdt(u0,t,tdir,dtmax,uEltype(abstol),uEltypeNoUnits(reltol),internalnorm,f!,order))
  end

  if sign(dt)!=tdir && dt!=tType(0)
    error("dt has the wrong sign. Exiting")
  end

  rate_prototype = u/zero(t)
  rateType = typeof(rate_prototype) ## Can be different if united

  saveat = tType[convert(tType,x) for x in setdiff(saveat,tspan)]

  ### Algorithm-specific defaults ###

  abstol = uEltype(1)*abstol

  isspecialdense(alg) ? ksEltype = Vector{rateType} : ksEltype = rateType

  # Have to convert incase passed in wrong.
  timeseries = convert(Vector{uType},timeseries_init)
  ts = convert(Vector{tType},ts_init)
  ks = convert(Vector{ksEltype},ks_init)

  copyat_or_push!(ts,1,t)
  copyat_or_push!(timeseries,1,u)

  if !isspecialdense(alg)
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

  calcprevs = calck || !(typeof(callback)<:Void) # Calculate the previous values
  tprev = t
  dtcache = tType(dt)
  dt_mod = tTypeNoUnits(1)
  iter = 0
  saveiter = 1 # Starts at 1 so first save is at 2
  saveiter_dense = 1
  cursaveat = 1
  kshortsize = 1
  reeval_fsal = false

  integrator = ODEIntegrator{algType,uType,tType,tTypeNoUnits,eltype(ks),typeof(sol),
                             typeof(rate_prototype),typeof(f!),
                             typeof(event_cache),typeof(opts)}(
                             sol,u,k,t,tType(dt),f!,uprev,kprev,tprev,
                             Ts,adaptiveorder,order,
                             alg,rate_prototype,
                             notsaveat_idxs,calcprevs,dtcache,dt_mod,tdir,
                             iter,saveiter,saveiter_dense,cursaveat,
                             event_cache,kshortsize,reeval_fsal,opts)
  integrator
end

function solve!(integrator::ODEIntegrator;timeseries_errors = true,dense_errors = false)
  @code_warntype ode_solve(integrator)
  ode_solve(integrator)

  if typeof(integrator.sol.prob) <: AbstractODETestProblem
    calculate_solution_errors!(integrator.sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
  end
  nothing
end
