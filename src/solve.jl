function solve(
  prob::AbstractODEProblem,
  alg::algType,timeseries=[],ts=[],ks=[],recompile::Type{Val{recompile_flag}}=Val{true};
  kwargs...) where {algType<:OrdinaryDiffEqAlgorithm,recompile_flag}

  integrator = init(prob,alg,timeseries,ts,ks,recompile;kwargs...)
  solve!(integrator)
  integrator.sol
end

function init(
  prob::AbstractODEProblem,
  alg::algType,timeseries_init=typeof(prob.u0)[],
  ts_init=eltype(prob.tspan)[],ks_init=[],
  recompile::Type{Val{recompile_flag}}=Val{true};
  timeseries_steps = 1,
  saveat = eltype(prob.tspan)[],
  tstops = eltype(prob.tspan)[],
  d_discontinuities= eltype(prob.tspan)[],
  save_idxs = nothing,
  save_everystep = isempty(saveat),
  save_timeseries = nothing,save_start = true,save_end = true,
  callback=nothing,
  dense = save_everystep && !(typeof(alg) <: FunctionMap),
  calck = (callback != nothing && callback != CallbackSet()) || # Empty callback
          (prob.callback != nothing && prob.callback != CallbackSet()) || # Empty prob.callback
          (!isempty(setdiff(saveat,tstops)) || dense), # and no dense output
  dt = typeof(alg) <: FunctionMap && isempty(tstops) ? eltype(prob.tspan)(1) : eltype(prob.tspan)(0),
  adaptive = isadaptive(alg),
  gamma=gamma_default(alg),
  abstol=nothing,
  reltol=nothing,
  qmax=qmax_default(alg),qmin=qmin_default(alg),
  qsteady_min = qsteady_min_default(alg),
  qsteady_max = qsteady_max_default(alg),
  qoldinit=1//10^4, fullnormalize=true,
  failfactor = 2,
  beta2=nothing,
  beta1=nothing,
  maxiters = 1000000,
  dtmax=eltype(prob.tspan)((prob.tspan[end]-prob.tspan[1])),
  dtmin= typeof(one(eltype(prob.tspan))) <: AbstractFloat ? 10*eps(eltype(prob.tspan)) :
         typeof(one(eltype(prob.tspan))) <: Integer ? 0 :
         eltype(prob.tspan)(1//10^(10)),
  internalnorm = ODE_DEFAULT_NORM,
  isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
  unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
  verbose = true, force_dtmin = false,
  timeseries_errors = true, dense_errors=false,
  advance_to_tstop = false,stop_at_next_tstop=false,
  initialize_save = true,
  progress=false,progress_steps=1000,progress_name="ODE",
  progress_message = ODE_DEFAULT_PROG_MESSAGE,
  userdata=nothing,
  allow_extrapolation = alg_extrapolates(alg),
  initialize_integrator=true,kwargs...) where {algType<:OrdinaryDiffEqAlgorithm,recompile_flag}

  if typeof(prob.f)<:Tuple
    if min((mm != I for mm in prob.mass_matrix)...)
      error("This solver is not able to use mass matrices.")
    end
  elseif !(typeof(prob)<:DiscreteProblem) &&
         !(typeof(alg) <:MassMatrixAlgorithms) &&
         prob.mass_matrix != I
    error("This solver is not able to use mass matrices.")
  end

  if !isempty(saveat) && dense
    warn("Dense output is incompatible with saveat. Please use the SavingCallback from the Callback Library to mix the two behaviors.")
  end

  if (eltype(prob.u0) <: Dual && !(eltype(prob.tspan)<:Dual) ||
     !(eltype(prob.u0) <: Dual) && eltype(prob.tspan)<:Dual) && adaptive
     warn("Autodifferentiation through the solver with adaptive timestepping requires both time and states to be dual numbers. Please see the FAQ.")
  end

  tType = eltype(prob.tspan)
  tspan = prob.tspan
  tdir = sign(tspan[end]-tspan[1])

  t = tspan[1]

  if (((!(typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm) && !(typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm)) || !adaptive) && dt == tType(0) && isempty(tstops)) && !(typeof(alg) <: FunctionMap)
      error("Fixed timestep methods require a choice of dt or choosing the tstops")
  end

  f = prob.f
  p = prob.p

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
  uBottomEltype = recursive_bottom_eltype(u)
  uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

  ks = Vector{uType}(0)

  uEltypeNoUnits = recursive_unitless_eltype(u)
  tTypeNoUnits   = typeof(one(tType))

  if typeof(alg) <: FunctionMap
    abstol_internal = real.(zero(u))
  elseif abstol == nothing
    if uBottomEltypeNoUnits == uBottomEltype || !(typeof(u) <: ArrayPartition)
      abstol_internal = real(uBottomEltype(uBottomEltype(1)*1//10^6))
    else
      abstol_internal = real.(ones(u).*1//10^6)
    end
  else
    abstol_internal = real.(abstol)
  end

  if typeof(alg) <: FunctionMap
    reltol_internal = real.(zero(first(u)/t))
  elseif reltol == nothing
    reltol_internal = real(uBottomEltypeNoUnits(1//10^3))
  else
    reltol_internal = real.(reltol)
  end

  dtmax > zero(dtmax) && tdir < 0 && (dtmax *= tdir) # Allow positive dtmax, but auto-convert
  # dtmin is all abs => does not care about sign already.

  if isinplace(prob) && typeof(u) <: AbstractArray && eltype(u) <: Number # Could this be more efficient for other arrays?
    if !(typeof(u) <: ArrayPartition)
      rate_prototype = similar(u,typeof(oneunit(uBottomEltype)/oneunit(tType)),indices(u))
    else
      rate_prototype = similar(u, typeof.(oneunit.(recursive_bottom_eltype.(u.x))./oneunit(tType))...)
    end
  else
    rate_prototype = u./oneunit(tType)
  end
  rateType = typeof(rate_prototype) ## Can be different if united

  tstops_internal, saveat_internal, d_discontinuities_internal =
    tstop_saveat_disc_handling(tstops,saveat,d_discontinuities,tdir,tspan,tType)

  callbacks_internal = CallbackSet(callback,prob.callback)


  ### Algorithm-specific defaults ###
  if save_idxs == nothing
    ksEltype = Vector{rateType}
  else
    ks_prototype = rate_prototype[save_idxs]
    ksEltype = Vector{typeof(ks_prototype)}
  end

  # Have to convert incase passed in wrong.
  if save_idxs == nothing
    timeseries = convert(Vector{uType},timeseries_init)
  else
    u_initial = u[save_idxs]
    timeseries = convert(Vector{typeof(u_initial)},timeseries_init)
  end
  ts = convert(Vector{tType},ts_init)
  ks = convert(Vector{ksEltype},ks_init)
  alg_choice = Int[]

  if !adaptive && save_everystep
    dt == 0 ? steps = length(tstops) : steps = round(Int,float((tspan[2]-tspan[1])/dt),RoundUp)
    sizehint!(timeseries,steps+1)
    sizehint!(ts,steps+1)
    sizehint!(ks,steps+1)
  elseif save_everystep
    sizehint!(timeseries,1000)
    sizehint!(ts,1000)
    sizehint!(ks,1000)
  elseif !isempty(saveat_internal)
    sizehint!(timeseries,length(saveat_internal)+1)
    sizehint!(ts,length(saveat_internal)+1)
    sizehint!(ks,length(saveat_internal)+1)
  else
    sizehint!(timeseries,2)
    sizehint!(ts,2)
    sizehint!(ks,2)
  end

  if save_start
    saveiter = 1 # Starts at 1 so first save is at 2
    saveiter_dense = 1
    copyat_or_push!(ts,1,t)
    if save_idxs == nothing
      copyat_or_push!(timeseries,1,u)
      copyat_or_push!(ks,1,[rate_prototype])
    else
      copyat_or_push!(timeseries,1,u_initial,Val{false})
      copyat_or_push!(ks,1,[ks_prototype])
    end
  else
    saveiter = 0 # Starts at 0 so first save is at 1
    saveiter_dense = 0
  end

  QT = tTypeNoUnits <: Integer ? typeof(qmin) : tTypeNoUnits

  progress ? (prog = Juno.ProgressBar(name=progress_name)) : prog = nothing

  k = rateType[]

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

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol_internal,p,calck,Val{isinplace(prob)})

  if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
    id = CompositeInterpolationData(f,timeseries,ts,ks,alg_choice,dense,cache)
    beta2 == nothing && ( beta2=beta2_default(alg.algs[cache.current]) )
    beta1 == nothing && ( beta1=beta1_default(alg.algs[cache.current],beta2) )
  else
    id = InterpolationData(f,timeseries,ts,ks,dense,cache)
    beta2 == nothing && ( beta2=beta2_default(alg) )
    beta1 == nothing && ( beta1=beta1_default(alg,beta2) )
  end

  opts = DEOptions{typeof(abstol_internal),typeof(reltol_internal),QT,tType,
                   typeof(internalnorm),typeof(callbacks_internal),typeof(isoutofdomain),
                   typeof(progress_message),typeof(unstable_check),typeof(tstops_internal),
                   typeof(d_discontinuities_internal),typeof(userdata),typeof(save_idxs),
                   typeof(maxiters),typeof(tstops),typeof(saveat),
                   typeof(d_discontinuities)}(
                       maxiters,timeseries_steps,save_everystep,adaptive,abstol_internal,
                       reltol_internal,QT(gamma),QT(qmax),
                       QT(qmin),QT(qsteady_max),
                       QT(qsteady_min),QT(failfactor),tType(dtmax),
                       tType(dtmin),internalnorm,save_idxs,tstops_internal,saveat_internal,
                       d_discontinuities_internal,
                       tstops,saveat,d_discontinuities,
                       userdata,progress,progress_steps,
                       progress_name,progress_message,timeseries_errors,dense_errors,
                       QT(beta1),QT(beta2),QT(qoldinit),dense,
                       save_start,save_end,callbacks_internal,isoutofdomain,
                       unstable_check,verbose,
                       calck,force_dtmin,advance_to_tstop,stop_at_next_tstop)

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
  kshortsize = 0
  reeval_fsal = false
  u_modified = false
  eigen_est = 1/oneunit(tType) # rate/state = (state/time)/state = 1/t units
  EEst = tTypeNoUnits(1)
  just_hit_tstop = false
  isout = false
  accept_step = false
  force_stepfail = false
  last_stepfail = false
  event_last_time = false
  dtchangeable = isdtchangeable(alg)
  q11 = tTypeNoUnits(1)
  success_iter = 0
  erracc = tTypeNoUnits(1)
  dtacc = tType(1)

  integrator = ODEIntegrator{algType,uType,tType,typeof(p),typeof(eigen_est),
                             QT,typeof(tdir),typeof(k),SolType,
                             FType,typeof(prog),cacheType,
                             typeof(opts),fsal_typeof(alg,rate_prototype)}(
                             sol,u,k,t,tType(dt),f,p,uprev,uprev2,tprev,
                             alg,dtcache,dtchangeable,
                             dtpropose,tdir,eigen_est,EEst,QT(qoldinit),q11,
                             erracc,dtacc,success_iter,
                             iter,saveiter,saveiter_dense,prog,cache,
                             kshortsize,force_stepfail,last_stepfail,
                             just_hit_tstop,event_last_time,accept_step,
                             isout,reeval_fsal,
                             u_modified,opts)
  if initialize_integrator
    initialize_callbacks!(integrator, initialize_save)
    initialize!(integrator,integrator.cache)
    save_start && typeof(alg) <: CompositeAlgorithm && copyat_or_push!(alg_choice,1,integrator.cache.current)
  end

  if integrator.dt == zero(integrator.dt) && integrator.opts.adaptive
    auto_dt_reset!(integrator)
    if sign(integrator.dt)!=integrator.tdir && integrator.dt!=tType(0) && !isnan(integrator.dt)
      error("Automatic dt setting has the wrong sign. Exiting. Please report this error.")
    end
    if isnan(integrator.dt)
      if verbose
        warn("Automatic dt set the starting dt as NaN, causing instability.")
      end
    end
  elseif integrator.opts.adaptive && integrator.dt > zero(integrator.dt) && integrator.tdir < 0
    integrator.dt *= integrator.tdir # Allow positive dt, but auto-convert
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
  integrator.sol = solution_new_retcode(integrator.sol,:Success)
  nothing
end

# Helpers

function tstop_saveat_disc_handling(tstops,saveat,d_discontinuities,tdir,tspan,tType)

  if isempty(d_discontinuities) && isempty(tstops) # TODO: Specialize more
    tstops_vec = [tspan[2]]
  else
    tstops_vec = vec(collect(tType,Iterators.filter(x->tdir*tspan[1]<tdir*x≤tdir*tspan[end],Iterators.flatten((tstops,d_discontinuities,tspan[end])))))
  end

  if tdir>0
    tstops_internal = binary_minheap(tstops_vec)
  else
    tstops_internal = binary_maxheap(tstops_vec)
  end

  if typeof(saveat) <: Number
    if (tspan[1]:saveat:tspan[end])[end] == tspan[end]
      saveat_vec = convert(Vector{tType},collect(tType,tspan[1]+saveat:saveat:tspan[end]))
    else
      saveat_vec = convert(Vector{tType},collect(tType,tspan[1]+saveat:saveat:(tspan[end]-saveat)))
    end
  elseif isempty(saveat)
    saveat_vec = saveat
  else
    saveat_vec = vec(collect(tType,Iterators.filter(x->tdir*tspan[1]<tdir*x<tdir*tspan[end],saveat)))
  end

  if tdir>0
    saveat_internal = binary_minheap(saveat_vec)
  else
    saveat_internal = binary_maxheap(saveat_vec)
  end

  d_discontinuities_vec = vec(collect(d_discontinuities))

  if tdir>0
    d_discontinuities_internal = binary_minheap(d_discontinuities_vec)
  else
    d_discontinuities_internal = binary_maxheap(d_discontinuities_vec)
  end
  tstops_internal,saveat_internal,d_discontinuities_internal
end

function initialize_callbacks!(integrator, initialize_save = true)
  t = integrator.t
  u = integrator.u
  callbacks = integrator.opts.callback
  integrator.u_modified = true

  u_modified = initialize!(callbacks,u,t,integrator)

  # if the user modifies u, we need to fix previous values before initializing
  # FSAL in order for the starting derivatives to be correct
  if u_modified

    if isinplace(integrator.sol.prob)
      recursivecopy!(integrator.uprev,integrator.u)
    else
      integrator.uprev = integrator.u
    end

    if alg_extrapolates(integrator.alg)
      if isinplace(integrator.sol.prob)
        recursivecopy!(integrator.uprev2,integrator.uprev)
      else
        integrator.uprev2 = integrator.uprev
      end
    end

    if initialize_save &&
      (any((c)->c.save_positions[2],callbacks.discrete_callbacks) ||
      any((c)->c.save_positions[2],callbacks.continuous_callbacks))
      savevalues!(integrator,true)
    end
  end

  # reset this as it is now handled so the integrators should proceed as normal
  integrator.u_modified = false
end
