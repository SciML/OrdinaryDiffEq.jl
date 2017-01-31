function solve{uType,tType,isinplace,algType<:OrdinaryDiffEqAlgorithm,F,recompile_flag}(
  prob::AbstractODEProblem{uType,tType,isinplace,F},
  alg::algType,timeseries=[],ts=[],ks=[],recompile::Type{Val{recompile_flag}}=Val{true};
  kwargs...)

  integrator = init(prob,alg,timeseries,ts,ks,recompile;kwargs...)
  solve!(integrator)
  integrator.sol
end

function init{uType,tType,isinplace,algType<:OrdinaryDiffEqAlgorithm,F,recompile_flag}(
  prob::AbstractODEProblem{uType,tType,isinplace,F},
  alg::algType,timeseries_init=uType[],ts_init=tType[],ks_init=[],
  recompile::Type{Val{recompile_flag}}=Val{true};
  dt = tType(0),save_timeseries = true,
  timeseries_steps = 1,
  dense = save_timeseries,
  saveat = tType[],tstops = tType[],d_discontinuities= tType[],
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
  unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
  verbose = true, force_dtmin = false,
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

  d_discontinuities_col = collect(d_discontinuities)

  if tdir>0
    tstops_internal = binary_minheap(convert(Vector{tType},append!(collect(tstops),d_discontinuities_col)))
  else
    tstops_internal = binary_maxheap(convert(Vector{tType},append!(collect(tstops),d_discontinuities_col)))
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

  uEltypeNoUnits = typeof(recursive_one(u))
  tTypeNoUnits   = typeof(recursive_one(t))

  if dt == zero(dt) && adaptive
    dt = tType(ode_determine_initdt(u0,t,tdir,dtmax,uEltype(abstol),uEltypeNoUnits(reltol),internalnorm,prob,order))
  end

  if sign(dt)!=tdir && dt!=tType(0)
    error("dt has the wrong sign. Exiting")
  end

  if typeof(u) <: AbstractArray
    rate_prototype = similar(u/zero(t),indices(u)) # rate doesn't need type info
  else
    rate_prototype = u/zero(t)
  end
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

  d_discontinuities_vec =  convert(Vector{tType},d_discontinuities_col)

  if tdir>0
    d_discontinuities_internal = binary_minheap(d_discontinuities_vec)
  else
    d_discontinuities_internal = binary_maxheap(d_discontinuities_vec)
  end

  callbacks_internal = CallbackSet(callback)


  ### Algorithm-specific defaults ###
  ksEltype = Vector{rateType}

  # Have to convert incase passed in wrong.
  timeseries = convert(Vector{uType},timeseries_init)
  ts = convert(Vector{tType},ts_init)
  ks = convert(Vector{ksEltype},ks_init)
  alg_choice = Int[]

  copyat_or_push!(ts,1,t)
  copyat_or_push!(timeseries,1,u)
  copyat_or_push!(ks,1,[rate_prototype])

  opts = DEOptions(Int(maxiters),timeseries_steps,save_timeseries,adaptive,uEltype(uEltype(1)*abstol),
    uEltypeNoUnits(reltol),tTypeNoUnits(gamma),tTypeNoUnits(qmax),tTypeNoUnits(qmin),
    dtmax,dtmin,internalnorm,
    tstops_internal,saveat_internal,d_discontinuities_internal,
    userdata,
    progress,progress_steps,
    progress_name,progress_message,
    timeseries_errors,dense_errors,
    tTypeNoUnits(beta1),tTypeNoUnits(beta2),tTypeNoUnits(qoldinit),dense,
    callbacks_internal,isoutofdomain,unstable_check,verbose,calck,force_dtmin,
    advance_to_tstop,stop_at_next_tstop)

  progress ? (prog = Juno.ProgressBar(name=progress_name)) : prog = nothing

  notsaveat_idxs = Int[1]

  k = ksEltype[]

  if uType <: Array
    uprev = copy(u)
  else
    uprev = deepcopy(u)
  end
  if alg_extrapolates(alg)
    if uType <: Array
      uprev2 = copy(u)
    else
      uprev2 = deepcopy(u)
    end
  else
    uprev2 = uprev
  end

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,Val{isinplace})

  if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
    id = CompositeInterpolationData(f,timeseries,ts,ks,alg_choice,notsaveat_idxs,cache)
  else
    id = InterpolationData(f,timeseries,ts,ks,notsaveat_idxs,cache)
  end

  if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
    sol = build_solution(prob,alg,ts,timeseries,
                      dense=dense,k=ks,interp=id,
                      alg_choice=alg_choice,
                      calculate_error = false)
  else
    sol = build_solution(prob,alg,ts,timeseries,
                      dense=dense,k=ks,interp=id,
                      calculate_error = false)
  end

  if recompile_flag == true
    FType = typeof(f)
    SolType = typeof(sol)
    cacheType = typeof(cache)
  else
    FType = Function
    SolType = AbstractODESolution
    cacheType =  OrdinaryDiffEqCache
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
  isout = false
  accept_step = false
  dtchangeable = isdtchangeable(alg)
  q11 = tTypeNoUnits(1)

  integrator = ODEIntegrator{algType,uType,tType,
                             tTypeNoUnits,typeof(tdir),eltype(ks),SolType,
                             typeof(rate_prototype),FType,typeof(prog),cacheType,
                             typeof(opts)}(
                             sol,u,k,t,tType(dt),f,uprev,uprev2,tprev,
                             alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,
                             dtpropose,dt_mod,tdir,EEst,qoldinit,q11,
                             iter,saveiter,saveiter_dense,prog,cache,
                             kshortsize,just_hit_tstop,accept_step,isout,reeval_fsal,u_modified,opts)
  if initialize_integrator
    initialize!(integrator,integrator.cache)
  end
  integrator
end

function solve!(integrator::ODEIntegrator)
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
