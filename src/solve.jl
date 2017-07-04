function solve{algType<:OrdinaryDiffEqAlgorithm,recompile_flag}(
  prob::AbstractODEProblem,
  alg::algType,timeseries=[],ts=[],ks=[],recompile::Type{Val{recompile_flag}}=Val{true};
  kwargs...)

  integrator = init(prob,alg,timeseries,ts,ks,recompile;kwargs...)
  solve!(integrator)
  integrator.sol
end

function init{algType<:OrdinaryDiffEqAlgorithm,recompile_flag}(
  prob::AbstractODEProblem,
  alg::algType,timeseries_init=typeof(prob.u0)[],ts_init=eltype(prob.tspan)[],ks_init=[],
  recompile::Type{Val{recompile_flag}}=Val{true};
  timeseries_steps = 1,
  saveat = eltype(prob.tspan)[],tstops = eltype(prob.tspan)[],d_discontinuities= eltype(prob.tspan)[],
  save_idxs = nothing,
  save_everystep = isempty(saveat),
  save_timeseries = nothing,save_start = true,
  dense = save_everystep && !(typeof(alg) <: Discrete),
  calck = (!isempty(setdiff(saveat,tstops)) || dense),
  dt = typeof(alg) <: Discrete && isempty(tstops) ? eltype(prob.tspan)(1) : eltype(prob.tspan)(0),
  adaptive = isadaptive(alg),
  gamma=9//10,
  abstol=nothing,
  reltol=nothing,
  qmax=qmax_default(alg),qmin=qmin_default(alg),
  qoldinit=1//10^4, fullnormalize=true,
  beta2=beta2_default(alg),
  beta1=beta1_default(alg,beta2),
  maxiters = 1000000,
  dtmax=eltype(prob.tspan)((prob.tspan[end]-prob.tspan[1])),
  dtmin=eltype(prob.tspan) <: AbstractFloat ? eltype(prob.tspan)(10)*eps(eltype(prob.tspan)) : eltype(prob.tspan)(1//10^(10)),
  internalnorm = ODE_DEFAULT_NORM,
  isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
  unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
  verbose = true, force_dtmin = false,
  timeseries_errors = true, dense_errors=false,
  advance_to_tstop = false,stop_at_next_tstop=false,
  progress=false,progress_steps=1000,progress_name="ODE",
  progress_message = ODE_DEFAULT_PROG_MESSAGE,
  userdata=nothing,callback=nothing,
  allow_extrapolation = alg_extrapolates(alg),
  initialize_integrator=true,kwargs...)

  if save_timeseries != nothing
    warn("save_timeseries is deprecated. Use save_everystep instead")
    save_everystep = save_timeseries
  end

  if typeof(prob.f)<:Tuple
    if min((mm != I for mm in prob.mass_matrix)...)
      error("This solver is not able to use mass matrices.")
    end
  elseif !(typeof(prob)<:DiscreteProblem) &&
         !(typeof(alg) <: Union{Rosenbrock23,Rosenbrock32}) &&
         prob.mass_matrix != I
    error("This solver is not able to use mass matrices.")
  end

  tType = eltype(prob.tspan)
  tspan = prob.tspan
  tdir = sign(tspan[end]-tspan[1])

  t = tspan[1]

  if ((!(typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm) && !(typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm)) || !adaptive) && dt == tType(0) && isempty(tstops)
      error("Fixed timestep methods require a choice of dt or choosing the tstops")
  end

  if tspan[1] == tspan[end]
    error("Timespan is trivial")
  end

  tstops_vec = collect(tType,Iterators.filter(x->tdir*tspan[1]<tdir*x≤tdir*tspan[end],Iterators.flatten((tstops,d_discontinuities,tspan[end]))))

  if tdir>0
    tstops_internal = binary_minheap(tstops_vec)
  else
    tstops_internal = binary_maxheap(tstops_vec)
  end

  f = prob.f

  # Get the control variables

  if typeof(prob.u0) <: Array
    u = recursivecopy(prob.u0)
  elseif typeof(prob.u0) <: Number
    u = prob.u0
  elseif typeof(prob.u0) <: Tuple
    u = ArrayPartition(prob.u0,Val{true})
  else
    u = deepcopy(prob.u0)
  end

  uType = typeof(u)
  uEltype = recursive_eltype(u)

  ks = Vector{uType}(0)

  order = alg_order(alg)

  uEltypeNoUnits = typeof(one(uEltype))
  tTypeNoUnits   = typeof(one(tType))

  if typeof(alg) <: Discrete
    abstol_internal = zero(u)
  elseif abstol == nothing
    if uEltypeNoUnits == uEltype || !(typeof(u) <: ArrayPartition)
      abstol_internal = uEltype(uEltype(1)*1//10^6)
    else
      abstol_internal = ones(u).*1//10^6
    end
  else
    abstol_internal = abstol
  end

  if typeof(alg) <: Discrete
    reltol_internal = zero(first(u)/t)
  elseif reltol == nothing
    reltol_internal = uEltypeNoUnits(1//10^3)
  else
    reltol_internal = reltol
  end

  dtmax > zero(dtmax) && tdir < 0 && (dtmax *= tdir) # Allow positive dtmax, but auto-convert
  # dtmin is all abs => does not care about sign already.
  if dt == zero(dt) && adaptive
    dt = tType(ode_determine_initdt(u,t,tdir,dtmax,abstol_internal,reltol_internal,internalnorm,prob,order))
    if sign(dt)!=tdir && dt!=tType(0)
      error("Automatic dt setting has the wrong sign. Exiting. Please report this error.")
    end
  elseif adaptive && dt > zero(dt) && tdir < 0
    dt *= tdir # Allow positive dt, but auto-convert
  end

  if isinplace(prob) && typeof(u) <: AbstractArray && eltype(u) <: Number # Could this be more efficient for other arrays?
    rate_prototype = similar(u,typeof(oneunit(uEltype)/oneunit(tType)))
  else
    rate_prototype = u./oneunit(tType)
  end
  rateType = typeof(rate_prototype) ## Can be different if united

  if typeof(saveat) <: Number
    saveat_vec = collect(tType,tspan[1]+saveat:saveat:(tspan[end]-saveat))
    # Exclude the endpoint because of floating point issues
  else
    saveat_vec = collect(tType,Iterators.filter(x->tdir*tspan[1]<tdir*x<tdir*tspan[end],saveat))
  end

  if tdir>0
    saveat_internal = binary_minheap(saveat_vec)
  else
    saveat_internal = binary_maxheap(saveat_vec)
  end

  d_discontinuities_vec =  collect(tType,d_discontinuities)

  if tdir>0
    d_discontinuities_internal = binary_minheap(d_discontinuities_vec)
  else
    d_discontinuities_internal = binary_maxheap(d_discontinuities_vec)
  end

  callbacks_internal = CallbackSet(callback,prob.callback)


  ### Algorithm-specific defaults ###
  ksEltype = Vector{rateType}

  # Have to convert incase passed in wrong.
  timeseries = convert(Vector{uType},timeseries_init)
  ts = convert(Vector{tType},ts_init)
  ks = convert(Vector{ksEltype},ks_init)
  alg_choice = Int[]

  if save_start
    saveiter = 1 # Starts at 1 so first save is at 2
    saveiter_dense = 1
    copyat_or_push!(ts,1,t)
    if save_idxs == nothing
      copyat_or_push!(timeseries,1,u)
    else
      copyat_or_push!(timeseries,1,u[save_idxs],Val{false})
    end
    copyat_or_push!(ks,1,[rate_prototype])
  else
    saveiter = 0 # Starts at 0 so first save is at 1
    saveiter_dense = 0
  end

  opts = DEOptions(maxiters,timeseries_steps,save_everystep,adaptive,abstol_internal,
    reltol_internal,tTypeNoUnits(gamma),tTypeNoUnits(qmax),tTypeNoUnits(qmin),
    tType(dtmax),tType(dtmin),internalnorm,save_idxs,
    tstops_internal,saveat_internal,d_discontinuities_internal,
    userdata,
    progress,progress_steps,
    progress_name,progress_message,
    timeseries_errors,dense_errors,
    tTypeNoUnits(beta1),tTypeNoUnits(beta2),tTypeNoUnits(qoldinit),dense,save_start,
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
  if allow_extrapolation
    if uType <: Array
      uprev2 = copy(u)
    else
      uprev2 = deepcopy(u)
    end
  else
    uprev2 = uprev
  end

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,Val{isinplace(prob)})

  if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
    id = CompositeInterpolationData(f,timeseries,ts,ks,alg_choice,notsaveat_idxs,dense,cache)
  else
    id = InterpolationData(f,timeseries,ts,ks,notsaveat_idxs,dense,cache)
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
  iter = 0
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
                             typeof(opts),fsal_typeof(alg,rate_prototype)}(
                             sol,u,k,t,tType(dt),f,uprev,uprev2,tprev,
                             alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,
                             dtpropose,tdir,EEst,qoldinit,q11,
                             iter,saveiter,saveiter_dense,prog,cache,
                             kshortsize,just_hit_tstop,accept_step,isout,reeval_fsal,u_modified,opts)
  if initialize_integrator
    initialize!(integrator,integrator.cache)
    initialize!(callbacks_internal,t,u,integrator)
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

  if typeof(integrator.sol.prob.f) <: Tuple
    f = integrator.sol.prob.f[1]
  else
    f = integrator.sol.prob.f
  end

  if has_analytic(f)
    calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
  end
  integrator.sol.retcode = :Success
  nothing
end
