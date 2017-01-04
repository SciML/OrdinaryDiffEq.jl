function solve{uType,tType,isinplace,algType<:OrdinaryDiffEqAlgorithm,F}(
  prob::AbstractODEProblem{uType,tType,isinplace,F},
  alg::algType,timeseries=[],ts=[],ks=[];
  kwargs...)

  integrator = init(prob,alg,timeseries,ts,ks;kwargs...)
  solve!(integrator)
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
  gamma=9//10,
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
  timeseries_errors = true, dense_errors=false,
  advance_to_tstop = false,stop_at_next_tstop=false,
  progress=false,progress_steps=1000,progress_name="ODE",
  progress_message = ODE_DEFAULT_PROG_MESSAGE,
  userdata=nothing,callback=nothing,kwargs...)

  tspan = prob.tspan
  tdir = sign(tspan[end]-tspan[1])

  t = tspan[1]

  if !(typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm) && dt == tType(0) && isempty(tstops)
      error("Fixed timestep methods require a choice of dt or choosing the tstops")
  end

  if !isempty(tstops) && tdir*tspan[end] < tdir*maximum(tstops)
      error("Final saving timepoint is past the solving timespan")
  end
  if !isempty(tstops) && tdir*t > tdir*minimum(tstops)
      error("First saving timepoint is before the solving timespan")
  end

  if tdir>0
    tstops_internal = binary_minheap(convert(Vector{tType},collect(tstops)))
  else
    tstops_internal = binary_maxheap(convert(Vector{tType},collect(tstops)))
  end

  if !isempty(tstops) && tstops[end] != tspan[2]
    push!(tstops_internal,tspan[2])
  elseif isempty(tstops)
    push!(tstops_internal,tspan[2])
  end

  if top(tstops_internal) == tspan[1]
    pop!(tstops_internal)
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

  saveat_vec =  convert(Vector{tType},collect(saveat))
  if !isempty(saveat_vec) && saveat_vec[end] == tspan[2]
    pop!(saveat_vec)
  end

  if tdir>0
    saveat_internal = binary_minheap(saveat_vec)
  else
    saveat_internal = binary_maxheap(saveat_vec)
  end

  if !isempty(saveat_internal) && top(saveat_internal) == tspan[1]
    pop!(saveat_internal)
  end





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
    uEltypeNoUnits(reltol),tTypeNoUnits(gamma),tTypeNoUnits(qmax),tTypeNoUnits(qmin),
    dtmax,dtmin,internalnorm,tstops_internal,saveat_internal,userdata,
    progress,progress_steps,
    progress_name,progress_message,
    timeseries_errors,dense_errors,
    tTypeNoUnits(beta1),tTypeNoUnits(beta2),tTypeNoUnits(qoldinit),dense,
    callback,isoutofdomain,calck,advance_to_tstop,stop_at_next_tstop)

  progress ? (prog = Juno.ProgressBar(name=progress_name)) : prog = nothing

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

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,uprev,kprev,f!,t)

  if dense
    id = InterpolationData(f!,timeseries,ts,ks,notsaveat_idxs)
    interp = (tvals) -> ode_interpolation(cache,tvals,id)
  else
    interp = (tvals) -> nothing
  end

  sol = build_solution(prob,alg,ts,timeseries,
                    dense=dense,k=ks,interp=interp,
                    calculate_error = false)

  calcprevs = calck || !(typeof(callback)<:Void) # Calculate the previous values
  tprev = t
  dtcache = tType(dt)
  dtpropose = tType(dt)
  dt_mod = tTypeNoUnits(1)
  iter = 0
  saveiter = 1 # Starts at 1 so first save is at 2
  saveiter_dense = 1
  kshortsize = 1
  reeval_fsal = false
  qminc = inv(qmin) #facc1
  qmaxc = inv(qmax) #facc2
  EEst = tTypeNoUnits(1)
  just_hit_tstop = false
  accept_step = false

  integrator = ODEIntegrator{algType,uType,tType,
                             tTypeNoUnits,eltype(ks),typeof(sol),
                             typeof(rate_prototype),typeof(f!),typeof(prog),typeof(cache),
                             typeof(opts)}(
                             sol,u,k,t,tType(dt),f!,uprev,kprev,tprev,
                             adaptiveorder,order,
                             alg,rate_prototype,notsaveat_idxs,calcprevs,dtcache,
                             dtpropose,dt_mod,tdir,qminc,qmaxc,EEst,qoldinit,
                             iter,saveiter,saveiter_dense,prog,cache,
                             kshortsize,just_hit_tstop,accept_step,reeval_fsal,opts)
  integrator
end

function solve!(integrator::ODEIntegrator)
  #@code_warntype ode_solve(integrator)
  #for i in integrator end
  ode_solve(integrator)

  if typeof(integrator.sol.prob) <: AbstractODETestProblem
    calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
  end
  nothing
end

function ode_solve(integrator::ODEIntegrator)
  initialize!(integrator,integrator.cache)
  @inbounds while !isempty(integrator.opts.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.opts.tstops)
      loopheader!(integrator)
      @ode_exit_conditions
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
      if isempty(integrator.opts.tstops)
        break
      end
    end
    !isempty(integrator.opts.tstops) && pop!(integrator.opts.tstops)
  end
  ode_postamble!(integrator)
  nothing
end
