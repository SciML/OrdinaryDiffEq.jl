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
  adaptive = isadaptive(alg),
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
  userdata=nothing,callback=nothing,
  initialize_integrator=true,kwargs...)

  tspan = prob.tspan
  tdir = sign(tspan[end]-tspan[1])

  t = tspan[1]

  if (!(typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm) && !(typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm)) && dt == tType(0) && isempty(tstops)
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
  f = prob.f
  u0 = prob.u0
  uEltype = eltype(u0)

  # Get the control variables

  (uType<:Array || uType <: Number) ? u = copy(u0) : u = deepcopy(u0)

  ks = Vector{uType}(0)

  order = alg_order(alg)

  if typeof(alg) <: ExplicitRK
    @unpack order = alg.tableau
  elseif (typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm) && typeof(alg.algs[1]) <: ExplicitRK
    @unpack order = alg.algs[1].tableau
  end

  uEltypeNoUnits = typeof(recursive_one(u))
  tTypeNoUnits   = typeof(recursive_one(t))

  if dt == zero(dt) && adaptive
    dt = tType(ode_determine_initdt(u0,t,tdir,dtmax,uEltype(abstol),uEltypeNoUnits(reltol),internalnorm,f,order))
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

  ksEltype = Vector{rateType}

  # Have to convert incase passed in wrong.
  timeseries = convert(Vector{uType},timeseries_init)
  ts = convert(Vector{tType},ts_init)
  ks = convert(Vector{ksEltype},ks_init)
  alg_choice = Int[]

  copyat_or_push!(ts,1,t)
  copyat_or_push!(timeseries,1,u)
  copyat_or_push!(ks,1,[rate_prototype])

  if typeof(callback) <: DECallback
    # Change it to a tuple
    callback_internal = (callback,)
  else
    callback_internal = callback
  end

  opts = DEOptions(maxiters,timeseries_steps,save_timeseries,adaptive,uEltype(abstol),
    uEltypeNoUnits(reltol),tTypeNoUnits(gamma),tTypeNoUnits(qmax),tTypeNoUnits(qmin),
    dtmax,dtmin,internalnorm,tstops_internal,saveat_internal,userdata,
    progress,progress_steps,
    progress_name,progress_message,
    timeseries_errors,dense_errors,
    tTypeNoUnits(beta1),tTypeNoUnits(beta2),tTypeNoUnits(qoldinit),dense,
    callback_internal,isoutofdomain,calck,advance_to_tstop,stop_at_next_tstop)

  progress ? (prog = Juno.ProgressBar(name=progress_name)) : prog = nothing

  notsaveat_idxs = Int[1]

  k = ksEltype[]

  if uType <: Array
    uprev = copy(u)
  else
    uprev = deepcopy(u)
  end

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,uprev,f,t,Val{isinplace})

  if dense
    if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
      id = CompositeInterpolationData(f,timeseries,ts,ks,alg_choice,notsaveat_idxs)
    else
      id = InterpolationData(f,timeseries,ts,ks,notsaveat_idxs)
    end
    interp = (tvals) -> ode_interpolation(cache,tvals,id)
  else
    interp = (tvals) -> nothing
  end

  if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
    sol = build_solution(prob,alg,ts,timeseries,
                      dense=dense,k=ks,interp=interp,
                      alg_choice=alg_choice,
                      calculate_error = false)
  else
    sol = build_solution(prob,alg,ts,timeseries,
                      dense=dense,k=ks,interp=interp,
                      calculate_error = false)
  end

  tprev = t
  dtcache = tType(dt)
  dtpropose = tType(dt)
  dt_mod = tTypeNoUnits(1)
  iter = 0
  saveiter = 1 # Starts at 1 so first save is at 2
  saveiter_dense = 1
  kshortsize = 1
  reeval_fsal = false
  u_modified = false
  EEst = tTypeNoUnits(1)
  just_hit_tstop = false
  accept_step = false
  dtchangeable = isdtchangeable(alg)
  q11 = tTypeNoUnits(1)

  integrator = ODEIntegrator{algType,uType,tType,
                             tTypeNoUnits,typeof(tdir),eltype(ks),typeof(sol),
                             typeof(rate_prototype),typeof(f),typeof(prog),typeof(cache),
                             typeof(opts)}(
                             sol,u,k,t,tType(dt),f,uprev,tprev,
                             alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,
                             dtpropose,dt_mod,tdir,EEst,qoldinit,q11,
                             iter,saveiter,saveiter_dense,prog,cache,
                             kshortsize,just_hit_tstop,accept_step,reeval_fsal,u_modified,opts)
  if initialize_integrator
    initialize!(integrator,integrator.cache)
  end
  integrator
end

function solve!(integrator::ODEIntegrator)
  # Should be:
  #for i in integrator end
  # But the performance is bad!
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
    handle_tstop!(integrator)
  end
  postamble!(integrator)
  if typeof(integrator.sol.prob) <: AbstractODETestProblem
    calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
  end
  nothing
end
