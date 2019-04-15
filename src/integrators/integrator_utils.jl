save_idxsinitialize(integrator,cache::OrdinaryDiffEqCache,::Type{uType}) where {uType} =
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
  choose_algorithm!(integrator,integrator.cache)
  fix_dt_at_bounds!(integrator)
  modify_dt_for_tstops!(integrator)
  integrator.force_stepfail = false
end

last_step_failed(integrator::ODEIntegrator) =
  integrator.last_stepfail && !integrator.opts.adaptive

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

function savevalues!(integrator::ODEIntegrator,force_save=false,reduce_size=true)::Tuple{Bool,Bool}
  saved, savedexactly = false, false
  !integrator.opts.save_on && return saved, savedexactly
  while !isempty(integrator.opts.saveat) && integrator.tdir*top(integrator.opts.saveat) <= integrator.tdir*integrator.t # Perform saveat
    integrator.saveiter += 1; saved = true
    curt = pop!(integrator.opts.saveat)
    if curt!=integrator.t # If <t, interpolate
      DiffEqBase.addsteps!(integrator)
      Θ = (curt - integrator.tprev)/integrator.dt
      val = ode_interpolant(Θ,integrator,integrator.opts.save_idxs,Val{0}) # out of place, but no force copy later
      copyat_or_push!(integrator.sol.t,integrator.saveiter,curt)
      save_val = val
      copyat_or_push!(integrator.sol.u,integrator.saveiter,save_val,Val{false})
      if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    else # ==t, just save
      savedexactly = true
      copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
      if integrator.opts.save_idxs === nothing
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
      else
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
      end
      if typeof(integrator.alg) <: FunctionMap || integrator.opts.dense
        integrator.saveiter_dense +=1
        if integrator.opts.dense
          if integrator.opts.save_idxs === nothing
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
  if force_save || integrator.opts.save_everystep
    integrator.saveiter += 1; saved, savedexactly = true, true
    if integrator.opts.save_idxs === nothing
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    else
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
    end
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if typeof(integrator.alg) <: FunctionMap || integrator.opts.dense
      integrator.saveiter_dense +=1
      if integrator.opts.dense
        if integrator.opts.save_idxs === nothing
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
  return saved, savedexactly
end

function postamble!(integrator::ODEIntegrator)
  solution_endpoint_match_cur_integrator!(integrator)
  resize!(integrator.sol.t,integrator.saveiter)
  resize!(integrator.sol.u,integrator.saveiter)
  resize!(integrator.sol.k,integrator.saveiter_dense)
  if integrator.opts.progress
    @logmsg(-1,
    integrator.opts.progress_name,
    _id = :OrdinaryDiffEq,
    message=integrator.opts.progress_message(integrator.dt,integrator.u,integrator.p,integrator.t),
    progress="done")
  end
end

function solution_endpoint_match_cur_integrator!(integrator)
  if integrator.opts.save_end && (integrator.saveiter == 0 || integrator.sol.t[integrator.saveiter] !=  integrator.t)
    integrator.saveiter += 1
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if integrator.opts.save_idxs === nothing
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    else
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
    end
    if typeof(integrator.alg) <: FunctionMap || integrator.opts.dense
      integrator.saveiter_dense +=1
      if integrator.opts.dense
        if integrator.opts.save_idxs === nothing
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
  if iszero(EEst)
    q = inv(integrator.opts.qmax)
  else
    @fastmath q11 = EEst^beta1
    @fastmath q = q11/(qold^beta2)
    integrator.q11 = q11
    @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),q/integrator.opts.gamma))
  end
  q
end

function step_accept_controller!(integrator,alg,q)
  if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    q = one(q)
  end
  integrator.qold = max(integrator.EEst,integrator.opts.qoldinit)
  integrator.dt/q #dtnew
end
function step_reject_controller!(integrator,alg)
  integrator.dt = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
end

const StandardControllerAlgs = Union{GenericImplicitEuler,GenericTrapezoid,VCABM}
#const NordAlgs = Union{AN5, JVODE}

function stepsize_controller!(integrator, alg::JVODE)
  #η = choose_η!(integrator, integrator.cache)
  if iszero(integrator.EEst)
    η = integrator.opts.qmax
  else
    η = integrator.cache.η
    integrator.qold = η
  end
  η
end
function step_accept_controller!(integrator,alg::JVODE,η)
  return η * integrator.dt  # dtnew
end
function step_reject_controller!(integrator,alg::JVODE)
  integrator.dt *= integrator.qold
end

function stepsize_controller!(integrator, alg::QNDF)
  cnt = integrator.iter
  if cnt <= 3
    # std controller
    if iszero(integrator.EEst)
      q = inv(integrator.opts.qmax)
    else
      qtmp = integrator.EEst^(1/(get_current_adaptive_order(integrator.alg,integrator.cache)+1))/integrator.opts.gamma
      @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp))
      integrator.qold = integrator.dt/q
    end
  else
    q = integrator.dt/integrator.cache.h
    integrator.qold = integrator.dt/q
  end
  q
end

function step_accept_controller!(integrator,alg::QNDF,q)
  return integrator.dt/q  # dtnew
end
function step_reject_controller!(integrator,alg::QNDF)
  integrator.dt = integrator.qold
end


function stepsize_controller!(integrator,alg::Union{StandardControllerAlgs,
                              OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,:Standard}}) where {CS, AD}
  # Standard stepsize controller
  if iszero(integrator.EEst)
    q = inv(integrator.opts.qmax)
  else
    qtmp = integrator.EEst^(1/(get_current_adaptive_order(integrator.alg,integrator.cache)+1))/integrator.opts.gamma
    @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp))
    integrator.qold = integrator.dt/q
  end
  q
end
function step_accept_controller!(integrator,alg::Union{StandardControllerAlgs,
                                 OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,:Standard}},q) where {CS, AD}
  integrator.dt/q # dtnew
end
function step_reject_controller!(integrator,alg::Union{StandardControllerAlgs,
                                 OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,:Standard}}) where {CS, AD}
  integrator.dt = integrator.qold
end

const PredictiveControllerAlgs = Union{RKC,RadauIIA5} # Union{RKC,OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,:Predictive}} where {CS, AD}
function stepsize_controller!(integrator,alg::PredictiveControllerAlgs)
  # Gustafsson predictive stepsize controller

  if iszero(integrator.EEst)
    q = inv(integrator.opts.qmax)
  else
    gamma = integrator.opts.gamma
    if typeof(alg) <: Union{RKC,IRKC}
      fac = gamma
    else
      if alg isa RadauIIA5
        @unpack nl_iters = integrator.cache
        @unpack max_iter = alg
      else
        @unpack nl_iters, max_iter = integrator.cache.nlsolve.cache
      end
      niters = nl_iters
      fac = min(gamma,(1+2*max_iter)*gamma/(niters+2*max_iter))
    end
    expo = 1/(get_current_adaptive_order(integrator.alg,integrator.cache)+1)
    qtmp = (integrator.EEst^expo)/fac
    @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp))
    integrator.qold = q
  end
  q
end
function step_accept_controller!(integrator,alg::PredictiveControllerAlgs,q)
  if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    q = one(q)
  end
  if integrator.success_iter > 0
    expo = 1/(get_current_adaptive_order(integrator.alg,integrator.cache)+1)
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
function step_reject_controller!(integrator,alg::PredictiveControllerAlgs)
  if integrator.success_iter == 0
    integrator.dt *= 0.1
  else
    integrator.dt = integrator.dt/integrator.qold
  end
end


@inline function stepsize_controller!(integrator,alg::ExtrapolationMidpointDeuflhard)
  # Dummy function
  # ExtrapolationMidpointDeuflhard's stepsize scaling is stored in the cache;
  # it is computed by  stepsize_controller_internal! (in perfom_step!) resp. stepsize_predictor!
  # (in step_accept_controller! and step_reject_controller!)
  zero(typeof(integrator.opts.qmax))
end

function stepsize_controller_internal!(integrator,alg::ExtrapolationMidpointDeuflhard)
  # Standard stepsize controller
  # Compute and save the stepsize scaling based on the latest error estimate of the current order
  if iszero(integrator.EEst)
    q = inv(integrator.opts.qmax)
  else
    # Update gamma and beta1
    integrator.opts.beta1 = typeof(integrator.opts.beta1)(1 // (2integrator.cache.n_curr + 1))
    integrator.opts.gamma = typeof(integrator.opts.gamma)(1 // 4)^integrator.opts.beta1
    # Compute new stepsize scaling
    qtmp = integrator.EEst^integrator.opts.beta1 / integrator.opts.gamma
    @fastmath q = max(inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), qtmp))
  end
  integrator.cache.Q[integrator.cache.n_curr - alg.n_min + 1] = q
end

function stepsize_predictor!(integrator,alg::ExtrapolationMidpointDeuflhard,n_new::Int64)
  # Compute and save the stepsize scaling for order n_new based on the latest error estimate of the current order.
  if iszero(integrator.EEst)
    q = inv(integrator.opts.qmax)
  else
    # Initialize
    @unpack t,EEst = integrator
    @unpack stage_number = integrator.cache
    tol = integrator.opts.internalnorm(integrator.opts.reltol,t) # Deuflhard's approach relies on EEstD ≈ ||relTol||
    s_curr = stage_number[integrator.cache.n_curr - alg.n_min + 1]
    s_new = stage_number[n_new - alg.n_min + 1]
    # Update gamma and beta1
    integrator.opts.beta1 = typeof(integrator.opts.beta1)(1 // (2integrator.cache.n_curr + 1))
    integrator.opts.gamma = typeof(integrator.opts.gamma)(1 // 4)^integrator.opts.beta1
    # Compute new stepsize scaling
    qtmp = (EEst * tol^(1.0 - s_curr / s_new))^integrator.opts.beta1 / integrator.opts.gamma
    @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp))
  end
  integrator.cache.Q[n_new - alg.n_min + 1] = q
end

function step_accept_controller!(integrator,alg::ExtrapolationMidpointDeuflhard,q)
  # Compute new order and stepsize, return new stepsize
  @unpack n_min, n_max = alg
  @unpack n_curr, n_old, Q = integrator.cache
  s = integrator.cache.stage_number

  # Compute new order based on available quantities
  tmp = (n_min:n_curr) .- n_min .+ 1 # Index range of quantities computed so far
  dt_new = Vector{eltype(Q)}(undef,length(tmp)+1)
  dt_new[1:end-1] = integrator.dt ./ Q[tmp] # Store for the possible new stepsizes
  dt_new[1:end-1] = max.(abs(integrator.opts.dtmin), min.(abs(integrator.opts.dtmax), abs.(dt_new[1:end-1]))) # Safety scaling

  # n_new is the most efficient order of the last step
  work = s[tmp] ./ dt_new[1:end-1]
  n_new = argmin(work) + n_min - 1

  # Check if n_new may be increased
  if n_new == n_curr < min(n_max, n_old + 1) # cf. win_max in perfom_step! of the last step
    # Predict stepsize scaling for order (n_new + 1)
    stepsize_predictor!(integrator, alg, n_new+1) # Update cache.Q

    # Compute and scale the corresponding stepsize
    dt_new[end] = integrator.dt ./ Q[tmp[end]+1]
    dt_new[end] = max(abs(integrator.opts.dtmin), min(abs(integrator.opts.dtmax), abs.(dt_new[end])))

    # Check if (n_new  + 1) would have been more efficient than n_new
    if work[end] > s[tmp[end]+1] / dt_new[end]
      n_new = n_new + 1
    end
  end

  integrator.cache.n_curr = n_new
  dt_new[n_new - n_min + 1]
end

function step_reject_controller!(integrator, alg::ExtrapolationMidpointDeuflhard)
  # Compute and save reduced stepsize dt_red of order n_old
  # Use the latest error estimate to predict dt_red if an estimate of order n_old is not available
  if integrator.cache.n_curr < integrator.cache.n_old
      stepsize_predictor!(integrator,alg,integrator.cache.n_old) # Update cache.Q
  end
  integrator.cache.n_curr = integrator.cache.n_old # Reset order for redoing the rejected step
  dt_red = integrator.dt / integrator.cache.Q[integrator.cache.n_old - integrator.alg.n_min + 1]
  dt_red = integrator.tdir*max(abs(integrator.opts.dtmin), min(abs(integrator.opts.dtmax), abs(dt_red))) # Safety scaling
  integrator.dt = dt_red
end

@inline function stepsize_controller!(integrator,alg::ExtrapolationMidpointHairerWanner)
  # Dummy function
  # ExtrapolationMidpointHairerWanner's stepsize scaling is stored in the cache;
  # it is computed by  stepsize_controller_internal! (in perfom_step!), step_accept_controller! or step_reject_controller!
  zero(typeof(integrator.opts.qmax))
end

function stepsize_controller_internal!(integrator,alg::ExtrapolationMidpointHairerWanner)
  # Standard stepsize controller
  # Compute and save the stepsize scaling based on the latest error estimate of the current order
  if iszero(integrator.EEst)
    q = inv(integrator.opts.qmax)
  else
    # Update gamma and beta1
    integrator.opts.beta1 = typeof(integrator.opts.beta1)(1 // (2integrator.cache.n_curr + 1))
    integrator.opts.gamma = typeof(integrator.opts.gamma)(65 // 100)^integrator.opts.beta1
    # Compute new stepsize scaling
    qtmp = integrator.EEst^integrator.opts.beta1 / integrator.opts.gamma
    @fastmath q = max(inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), qtmp))
  end
  integrator.cache.Q[integrator.cache.n_curr + 1] = q
end

function step_accept_controller!(integrator,alg::ExtrapolationMidpointHairerWanner,q)
  # Compute new order and stepsize, return new stepsize
  @unpack n_min, n_max = alg
  @unpack n_curr, n_old, Q, sigma = integrator.cache
  s = integrator.cache.stage_number

  # Compute new order based on available quantities
  win_min_old = n_old - 1 # cf. win_min in perfom_step! of the last step
  tmp = win_min_old:(n_curr + 1) # Index range for the new order
  dt_new = fill(zero(eltype(Q)),n_max+1)
  dt_new[tmp] = integrator.dt ./ Q[tmp] # Store for the possible new stepsizes
  dt_new[tmp] = max.(abs(integrator.opts.dtmin), min.(abs(integrator.opts.dtmax), abs.(dt_new[tmp]))) # Safety scaling
  work= Vector{eltype(Q)}(undef,n_max+1) # work[n] is the work for order (n-1)
  work[tmp] = s[tmp] ./ dt_new[tmp]
  # Order selection
  n_new = n_old
  if n_curr == n_min # Enforce n_min + 1 ≦ n_new
    n_new = n_min + 1
  else
    if n_curr <= n_old
      if work[n_curr-1] < sigma * work[n_curr]
        n_new = max(n_curr-1,n_old-1,n_min+1) # Enforce n_min + 1 ≦ n_new
      elseif work[n_curr] < sigma * work[n_curr-1]
        n_new = min(n_curr+1,n_max-1) # Enforce n_new ≦ n_max - 1
      else
        n_new = n_curr # n_min + 1 ≦ n_curr
      end
    else
      if work[n_old] < sigma *  work[n_old+1]
        n_new = max(n_old-1,n_min+1)  # Enforce n_min + 1 ≦ n_new
      end
      if work[n_curr+1] <  sigma * work[n_new+1]
        n_new = min(n_new+1,n_max-1) # Enforce n_new ≦ n_max - 1
      end
    end
  end
  integrator.cache.n_curr = n_new

  # Stepsize selection
  if n_new == n_curr + 1
    # Compute the new stepsize of order n_new based on the optimal stepsize of order n_curr
    dt_new[n_new+1] = s[n_curr + 2]/s[n_curr + 1 ] * dt_new[n_curr+1]
    dt_new[n_new+1] = max(abs(integrator.opts.dtmin), min(abs(integrator.opts.dtmax), abs(dt_new[n_new+1])))
  end
  dt_new[n_new + 1]
end

function step_reject_controller!(integrator, alg::ExtrapolationMidpointHairerWanner)
  # Compute and save order and stepsize for redoing the current step
  @unpack n_old, n_curr, Q = integrator.cache

  # Order selection
  n_red = n_old
  if n_curr == n_old - 1
    n_red = max(alg.n_min+1,n_old-1) # Enforce n_min + 1 ≦ n_red
  end
  integrator.cache.n_curr = n_red

  # Stepsize selection
  dt_red = integrator.dt / Q[n_red + 1]
  dt_red = integrator.tdir*max(abs(integrator.opts.dtmin), min(abs(integrator.opts.dtmax), abs(dt_red))) # Safety scaling
  integrator.dt = dt_red
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
      integrator.destats.naccept += 1
      integrator.last_stepfail = false
      dtnew = step_accept_controller!(integrator,integrator.alg,q)
      integrator.tprev = integrator.t
      # integrator.EEst has unitless type of integrator.t
      if typeof(integrator.EEst)<: AbstractFloat && !isempty(integrator.opts.tstops)
        tstop = top(integrator.opts.tstops)
        abs(ttmp - tstop) < 10eps(max(integrator.t,tstop)/oneunit(integrator.t))*oneunit(integrator.t) ?
                                  (integrator.t = tstop) : (integrator.t = ttmp)
      else
        integrator.t = ttmp
      end
      calc_dt_propose!(integrator,dtnew)
      handle_callbacks!(integrator)
    else # Reject
      integrator.destats.nreject += 1
    end
  elseif !integrator.opts.adaptive #Not adaptive
    integrator.destats.naccept += 1
    integrator.tprev = integrator.t
    # integrator.EEst has unitless type of integrator.t
    if typeof(integrator.EEst)<: AbstractFloat && !isempty(integrator.opts.tstops)
      tstop = top(integrator.opts.tstops)
      abs(ttmp - tstop) < 10eps(integrator.t/oneunit(integrator.t))*oneunit(integrator.t) ?
                                  (integrator.t = tstop) : (integrator.t = ttmp)
    else
      integrator.t = ttmp
    end
    integrator.last_stepfail = false
    integrator.accept_step = true
    integrator.dtpropose = integrator.dt
    handle_callbacks!(integrator)
  end
  if integrator.opts.progress && integrator.iter%integrator.opts.progress_steps==0
    @logmsg(-1,
    integrator.opts.progress_name,
    _id = :OrdinaryDiffEq,
    message=integrator.opts.progress_message(integrator.dt,integrator.u,integrator.p,integrator.t),
    progress=integrator.t/integrator.sol.prob.tspan[2])
  end
  (integrator.cache isa CompositeCache && integrator.eigen_est > integrator.destats.maxeig) && (integrator.destats.maxeig = integrator.eigen_est)
  nothing
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
              DiffEqBase.find_first_continuous_callback(integrator,continuous_callbacks...)
    if event_occurred
      integrator.event_last_time = idx
      continuous_modified,saved_in_cb = DiffEqBase.apply_callback!(integrator,continuous_callbacks[idx],time,upcrossing)
    else
      integrator.event_last_time = 0
    end
  end
  if !integrator.force_stepfail && !(typeof(discrete_callbacks)<:Tuple{})
    discrete_modified,saved_in_cb = DiffEqBase.apply_discrete_callback!(integrator,discrete_callbacks...)
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
      get_current_isfsal(integrator.alg, integrator.cache) && reset_fsal!(integrator)
  elseif get_current_isfsal(integrator.alg, integrator.cache)
    if integrator.reeval_fsal || integrator.u_modified || (typeof(integrator.alg)<:DP8 && !integrator.opts.calck) || (typeof(integrator.alg)<:Union{Rosenbrock23,Rosenbrock32} && !integrator.opts.adaptive)
        reset_fsal!(integrator)
    else # Do not reeval_fsal, instead copyto! over
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
        DiffEqBase.change_t_via_interpolation!(integrator, pop!(tstops), Val{true})
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

nlsolve_f(f, alg) = f isa SplitFunction && alg isa SplitAlgorithms ? f.f1 : f
nlsolve_f(integrator) =
  nlsolve_f(integrator.f, unwrap_alg(integrator, true))

function (integrator::ODEIntegrator)(t,deriv::Type=Val{0};idxs=nothing)
  current_interpolant(t,integrator,idxs,deriv)
end

(integrator::ODEIntegrator)(val::AbstractArray,t::Union{Number,AbstractArray},deriv::Type=Val{0};idxs=nothing) = current_interpolant!(val,t,integrator,idxs,deriv)
