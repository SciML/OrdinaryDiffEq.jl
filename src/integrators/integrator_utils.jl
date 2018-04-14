save_idxsinitialize{uType}(integrator,cache::OrdinaryDiffEqCache,::Type{uType}) =
                error("This algorithm does not have an initialization function")

function loopheader!(integrator)
  # Apply right after iterators / callbacks

  # Accept or reject the step
  if integrator.iter > 0
    if ((integrator.opts.adaptive && integrator.accept_step) || !integrator.opts.adaptive) && !integrator.force_stepfail
      integrator.success_iter += 1
      apply_step!(integrator)
    elseif integrator.opts.adaptive && !integrator.accept_step
      if integrator.isout
        integrator.dt = integrator.dt*integrator.opts.qmin
      elseif !integrator.force_stepfail
        step_reject_controller!(integrator,integrator.alg)
      end
    end
  end

  integrator.iter += 1
  fix_dt_at_bounds!(integrator)
  modify_dt_for_tstops!(integrator)
  integrator.force_stepfail = false
  choose_algorithm!(integrator,integrator.cache)
end

last_step_failed(integrator::ODEIntegrator) =
  integrator.last_stepfail && !integrator.opts.adaptive

@def ode_exit_conditions begin
  if check_error!(integrator) != :Success
    return integrator.sol
  end
end

function modify_dt_for_tstops!(integrator)
  tstops = integrator.opts.tstops
  if !isempty(tstops)
    if integrator.opts.adaptive
      if integrator.tdir > 0
        integrator.dt = min(abs(integrator.dt),abs(top(tstops)-integrator.t)) # step! to the end
      else
        integrator.dt = -min(abs(integrator.dt),abs(top(tstops)-integrator.t))
      end
    elseif integrator.dtcache == zero(integrator.t) && integrator.dtchangeable
      # Use integrator.opts.tstops
      integrator.dt = integrator.tdir*abs(top(tstops)-integrator.t)
  elseif integrator.dtchangeable && !integrator.force_stepfail
      # always try to step! with dtcache, but lower if a tstops
      # however, if force_stepfail then don't set to dtcache, and no tstop worry
      integrator.dt = integrator.tdir*min(abs(integrator.dtcache),abs(top(tstops)-integrator.t)) # step! to the end
    end
  end
end

function savevalues!(integrator::ODEIntegrator,force_save=false,reduce_size=true)
  while !isempty(integrator.opts.saveat) && integrator.tdir*top(integrator.opts.saveat) <= integrator.tdir*integrator.t # Perform saveat
    integrator.saveiter += 1
    curt = pop!(integrator.opts.saveat)
    if curt!=integrator.t # If <t, interpolate
      ode_addsteps!(integrator)
      Θ = (curt - integrator.tprev)/integrator.dt
      val = ode_interpolant(Θ,integrator,integrator.opts.save_idxs,Val{0}) # out of place, but no force copy later
      copyat_or_push!(integrator.sol.t,integrator.saveiter,curt)
      save_val = val
      copyat_or_push!(integrator.sol.u,integrator.saveiter,save_val,Val{false})
      if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    else # ==t, just save
      copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
      if integrator.opts.save_idxs ==nothing
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
      else
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
      end
      if typeof(integrator.alg) <: FunctionMap || integrator.opts.dense
        integrator.saveiter_dense +=1
        if integrator.opts.dense
          if integrator.opts.save_idxs ==nothing
            copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
          else
            copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,[k[integrator.opts.save_idxs] for k in integrator.k],Val{false})
          end
        end
      end
      if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    end
  end
  if force_save || (integrator.opts.save_everystep && integrator.iter%integrator.opts.timeseries_steps==0)
    integrator.saveiter += 1
    if integrator.opts.save_idxs == nothing
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    else
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
    end
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if typeof(integrator.alg) <: FunctionMap || integrator.opts.dense
      integrator.saveiter_dense +=1
      if integrator.opts.dense
        if integrator.opts.save_idxs == nothing
          copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
        else
          copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,[k[integrator.opts.save_idxs] for k in integrator.k],Val{false})
        end
      end
    end
    if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
      copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
    end
  end
  reduce_size && resize!(integrator.k,integrator.kshortsize)
end

function postamble!(integrator::ODEIntegrator)
  solution_endpoint_match_cur_integrator!(integrator)
  resize!(integrator.sol.t,integrator.saveiter)
  resize!(integrator.sol.u,integrator.saveiter)
  resize!(integrator.sol.k,integrator.saveiter_dense)
  !(typeof(integrator.prog)<:Void) && Juno.done(integrator.prog)
end

function solution_endpoint_match_cur_integrator!(integrator)
  if integrator.opts.save_end && (integrator.saveiter == 0 || integrator.sol.t[integrator.saveiter] !=  integrator.t)
    integrator.saveiter += 1
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if integrator.opts.save_idxs == nothing
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    else
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
    end
    if typeof(integrator.alg) <: FunctionMap || integrator.opts.dense
      integrator.saveiter_dense +=1
      if integrator.opts.dense
        if integrator.opts.save_idxs == nothing
          copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
        else
          copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,[k[integrator.opts.save_idxs] for k in integrator.k],Val{false})
        end
      end
    end
    if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
      copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
    end
  end
end

### Default is PI-controller
function stepsize_controller!(integrator,alg)
  # PI-controller
  EEst,beta1,q11,qold,beta2 = integrator.EEst, integrator.opts.beta1, integrator.q11,integrator.qold,integrator.opts.beta2
  @fastmath q11 = EEst^beta1
  @fastmath q = q11/(qold^beta2)
  integrator.q11 = q11
  @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),q/integrator.opts.gamma))
  if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    q = one(q)
  end
  q
end
function step_accept_controller!(integrator,alg,q)
  integrator.qold = max(integrator.EEst,integrator.opts.qoldinit)
  integrator.dt/q #dtnew
end
function step_reject_controller!(integrator,alg)
  integrator.dt = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
end

const StandardControllerAlgs = Union{GenericImplicitEuler,GenericTrapezoid}

function stepsize_controller!(integrator,alg::Union{StandardControllerAlgs,
                              OrdinaryDiffEqNewtonAdaptiveAlgorithm{:Standard}})
  # Standard stepsize controller
  qtmp = integrator.EEst^(1/(alg_adaptive_order(integrator.alg)+1))/integrator.opts.gamma
  @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp))
  integrator.qold = integrator.dt/q
  q
end
function step_accept_controller!(integrator,alg::Union{StandardControllerAlgs,
                              OrdinaryDiffEqNewtonAdaptiveAlgorithm{:Standard}},q)
  integrator.dt/q # dtnew
end
function step_reject_controller!(integrator,alg::Union{StandardControllerAlgs,
                              OrdinaryDiffEqNewtonAdaptiveAlgorithm{:Standard}})
  integrator.dt = integrator.qold
end

function stepsize_controller!(integrator,
                        alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm{:Predictive})

  # Gustafsson predictive stepsize controller
  gamma = integrator.opts.gamma
  niters = integrator.cache.newton_iters
  fac = min(gamma,(1+2*integrator.alg.max_newton_iter)*gamma/(niters+2*integrator.alg.max_newton_iter))
  expo = 1/(alg_order(integrator.alg)+1)
  qtmp = (integrator.EEst^expo)/fac
  @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp))
  if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    q = one(q)
  end
  integrator.qold = q
  q
end
function step_accept_controller!(integrator,
                      alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm{:Predictive},q)
  if integrator.success_iter > 0
    expo = 1/(alg_adaptive_order(integrator.alg)+1)
    qgus=(integrator.dtacc/integrator.dt)*(((integrator.EEst^2)/integrator.erracc)^expo)
    qgus = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qgus/integrator.opts.gamma))
    qacc=max(q,qgus)
  else
    qacc = q
  end
  integrator.dtacc = integrator.dt
  integrator.erracc = max(1e-2,integrator.EEst)
  integrator.dt/qacc
end
function step_reject_controller!(integrator,
                        alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm{:Predictive})
  if integrator.success_iter == 0
    integrator.dt *= 0.1
  else
    integrator.dt = integrator.dt/integrator.qold
  end
end

function loopfooter!(integrator)

  # Carry-over from callback
  # This is set to true if u_modified requires callback FSAL reset
  # But not set to false when reset so algorithms can check if reset occurred
  integrator.reeval_fsal = false
  integrator.u_modified = false

  ttmp = integrator.t + integrator.dt
  if integrator.force_stepfail
      if integrator.opts.adaptive
        integrator.dt = integrator.dt/integrator.opts.failfactor
      elseif integrator.last_stepfail
        return
      end
      integrator.last_stepfail = true
      integrator.accept_step = false
  elseif integrator.opts.adaptive
    q = stepsize_controller!(integrator,integrator.alg)
    integrator.isout = integrator.opts.isoutofdomain(integrator.u,integrator.p,ttmp)
    integrator.accept_step = (!integrator.isout && integrator.EEst <= 1.0) || (integrator.opts.force_dtmin && abs(integrator.dt) <= abs(integrator.opts.dtmin))
    if integrator.accept_step # Accept
      integrator.last_stepfail = false
      dtnew = step_accept_controller!(integrator,integrator.alg,q)
      integrator.tprev = integrator.t
      # integrator.EEst has unitless type of integrator.t
      if typeof(integrator.EEst)<: AbstractFloat && !isempty(integrator.opts.tstops)
        tstop = top(integrator.opts.tstops)
        abs(ttmp - tstop) < 10eps(typeof(integrator.t))*oneunit(integrator.t) ?
                                  (integrator.t = tstop) : (integrator.t = ttmp)
      else
        integrator.t = ttmp
      end
      calc_dt_propose!(integrator,dtnew)
      handle_callbacks!(integrator)
    end
  elseif !integrator.opts.adaptive #Not adaptive
    integrator.tprev = integrator.t
    # integrator.EEst has unitless type of integrator.t
    if typeof(integrator.EEst)<: AbstractFloat && !isempty(integrator.opts.tstops)
      tstop = top(integrator.opts.tstops)
      abs(ttmp - tstop) < 10eps(typeof(integrator.t))*oneunit(integrator.t) ?
                                  (integrator.t = tstop) : (integrator.t = ttmp)
    else
      integrator.t = ttmp
    end
    integrator.last_stepfail = false
    integrator.accept_step = true
    integrator.dtpropose = integrator.dt
    handle_callbacks!(integrator)
  end
  if !(typeof(integrator.prog)<:Void) && integrator.opts.progress && integrator.iter%integrator.opts.progress_steps==0
    Juno.msg(integrator.prog,integrator.opts.progress_message(integrator.dt,integrator.u,integrator.p,integrator.t))
    Juno.progress(integrator.prog,integrator.t/integrator.sol.prob.tspan[2])
  end
end

function handle_callbacks!(integrator)
  discrete_callbacks = integrator.opts.callback.discrete_callbacks
  continuous_callbacks = integrator.opts.callback.continuous_callbacks
  atleast_one_callback = false

  continuous_modified = false
  discrete_modified = false
  saved_in_cb = false
  if !(typeof(continuous_callbacks)<:Tuple{})
    time,upcrossing,event_occurred,idx,counter =
              find_first_continuous_callback(integrator,continuous_callbacks...)
    if event_occurred
      integrator.event_last_time = true
      continuous_modified,saved_in_cb = apply_callback!(integrator,continuous_callbacks[idx],time,upcrossing)
    else
      integrator.event_last_time = false
    end
  end
  if !integrator.force_stepfail && !(typeof(discrete_callbacks)<:Tuple{})
    discrete_modified,saved_in_cb = apply_discrete_callback!(integrator,discrete_callbacks...)
  end
  if !saved_in_cb
    savevalues!(integrator)
  end

  integrator.u_modified = continuous_modified || discrete_modified
  if integrator.u_modified
    handle_callback_modifiers!(integrator)
  end
end

function handle_callback_modifiers!(integrator::ODEIntegrator)
  integrator.reeval_fsal = true
end

function apply_step!(integrator)

  integrator.accept_step = false # yay we got here, don't need this no more

  #Update uprev
  if alg_extrapolates(integrator.alg)
    if isinplace(integrator.sol.prob)
      recursivecopy!(integrator.uprev2,integrator.uprev)
    else
      integrator.uprev2 = integrator.uprev
    end
  end
  if isinplace(integrator.sol.prob)
    recursivecopy!(integrator.uprev,integrator.u)
  else
    integrator.uprev = integrator.u
  end

  #Update dt if adaptive or if fixed and the dt is allowed to change
  if integrator.opts.adaptive || integrator.dtchangeable
    integrator.dt = integrator.dtpropose
  elseif integrator.dt != integrator.dtpropose && !integrator.dtchangeable
    error("The current setup does not allow for changing dt.")
  end

  # Update fsal if needed
  if !isempty(integrator.opts.d_discontinuities) &&
      top(integrator.opts.d_discontinuities) == integrator.t

      handle_discontinuities!(integrator)
      isfsal(integrator.alg) && reset_fsal!(integrator)
  elseif isfsal(integrator.alg)
    if integrator.reeval_fsal || integrator.u_modified || (typeof(integrator.alg)<:DP8 && !integrator.opts.calck) || (typeof(integrator.alg)<:Union{Rosenbrock23,Rosenbrock32} && !integrator.opts.adaptive)
        reset_fsal!(integrator)
    else # Do not reeval_fsal, instead copy! over
      if isinplace(integrator.sol.prob)
        recursivecopy!(integrator.fsalfirst,integrator.fsallast)
      else
        integrator.fsalfirst = integrator.fsallast
      end
    end
  end
end

function handle_discontinuities!(integrator)
    pop!(integrator.opts.d_discontinuities)
end

function calc_dt_propose!(integrator,dtnew)
  dtpropose = integrator.tdir*min(abs(integrator.opts.dtmax),abs(dtnew))
  dtpropose = integrator.tdir*max(abs(dtpropose),abs(integrator.opts.dtmin))
  integrator.dtpropose = dtpropose
end

function fix_dt_at_bounds!(integrator)
  if integrator.tdir > 0
    integrator.dt = min(integrator.opts.dtmax,integrator.dt)
  else
    integrator.dt = max(integrator.opts.dtmax,integrator.dt)
  end
  if integrator.tdir > 0
    integrator.dt = max(integrator.dt,integrator.opts.dtmin) #abs to fix complex sqrt issue at end
  else
    integrator.dt = min(integrator.dt,integrator.opts.dtmin) #abs to fix complex sqrt issue at end
  end
end

function handle_tstop!(integrator)
  tstops = integrator.opts.tstops
  if !isempty(tstops)
    t = integrator.t
    ts_top = top(tstops)
    if t == ts_top
      pop!(tstops)
      integrator.just_hit_tstop = true
    elseif integrator.tdir*t > integrator.tdir*ts_top
      if !integrator.dtchangeable
        change_t_via_interpolation!(integrator, pop!(tstops), Val{true})
        integrator.just_hit_tstop = true
      else
        error("Something went wrong. Integrator stepped past tstops but the algorithm was dtchangeable. Please report this error.")
      end
    end
  end
end

function reset_fsal!(integrator)
  # Under these condtions, these algorithms are not FSAL anymore
  if typeof(integrator.cache) <: OrdinaryDiffEqMutableCache ||
     (typeof(integrator.cache) <: CompositeCache &&
      typeof(integrator.cache.caches[1]) <: OrdinaryDiffEqMutableCache)
    integrator.f(integrator.fsalfirst,integrator.u,integrator.p,integrator.t)
  else
    integrator.fsalfirst = integrator.f(integrator.u,integrator.p,integrator.t)
  end
  # Do not set false here so it can be checked in the algorithm
  # integrator.reeval_fsal = false
end

function (integrator::ODEIntegrator)(t,deriv::Type=Val{0};idxs=nothing)
  current_interpolant(t,integrator,idxs,deriv)
end

(integrator::ODEIntegrator)(val::AbstractArray,t::Union{Number,AbstractArray},deriv::Type=Val{0};idxs=nothing) = current_interpolant!(val,t,integrator,idxs,deriv)
