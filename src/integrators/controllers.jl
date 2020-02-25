@inline function predictive_stepsize_controller!(integrator, alg)
  # Gustafsson predictive stepsize controller
  if iszero(integrator.EEst)
    q = inv(integrator.opts.qmax)
  else
    gamma = integrator.opts.gamma
    if typeof(alg) <: Union{RKC,IRKC,SERK2}
      fac = gamma
    else
      if alg isa RadauIIA5
        @unpack iter = integrator.cache
        @unpack maxiters = alg
      else
        @unpack iter, maxiters = integrator.cache.nlsolver
      end
      fac = min(gamma,(1+2*maxiters)*gamma/(iter+2*maxiters))
    end
    expo = 1/(get_current_adaptive_order(integrator.alg,integrator.cache)+1)
    qtmp = DiffEqBase.fastpow(integrator.EEst,expo)/fac
    @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp))
    integrator.qold = q
  end
  q
end

@inline function standard_stepsize_controller!(integrator, alg)
  # Standard stepsize controller
  if iszero(integrator.EEst)
    q = inv(integrator.opts.qmax)
  else
    qtmp = DiffEqBase.fastpow(integrator.EEst,1/(get_current_adaptive_order(integrator.alg,integrator.cache)+1))/integrator.opts.gamma
    @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp))
    integrator.qold = integrator.dt/q
  end
  q
end

@inline function PI_stepsize_controller!(integrator, alg)
  # PI-controller
  EEst,beta1,q11,qold,beta2 = integrator.EEst, integrator.opts.beta1, integrator.q11,integrator.qold,integrator.opts.beta2
  if iszero(EEst)
    q = inv(integrator.opts.qmax)
  else
    q11 = DiffEqBase.fastpow(EEst,beta1)
    q = q11/DiffEqBase.fastpow(qold,beta2)
    integrator.q11 = q11
    @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),q/integrator.opts.gamma))
  end
  q
end

@inline function stepsize_controller!(integrator,alg)
  if ispredictive(alg)
    predictive_stepsize_controller!(integrator, alg)
  elseif isstandard(alg)
    standard_stepsize_controller!(integrator, alg)
  else # Default is PI-controller
    PI_stepsize_controller!(integrator, alg)
  end
end

@inline function predictive_step_accept_controller!(integrator, alg, q)
  if integrator.success_iter > 0
    expo = 1/(get_current_adaptive_order(integrator.alg,integrator.cache)+1)
    qgus=(integrator.dtacc/integrator.dt)*DiffEqBase.fastpow((integrator.EEst^2)/integrator.erracc,expo)
    qgus = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qgus/integrator.opts.gamma))
    qacc=max(q,qgus)
  else
    qacc = q
  end
  if integrator.opts.qsteady_min <= qacc <= integrator.opts.qsteady_max
    qacc = one(qacc)
  end
  integrator.dtacc = integrator.dt
  integrator.erracc = max(1e-2,integrator.EEst)
  return integrator.dt/qacc
end

@inline function standard_step_accept_controller!(integrator, alg, q)
  if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    q = one(q)
  end
  integrator.dt/q
end

@inline function PI_step_accept_controller!(integrator, alg, q)
  if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    q = one(q)
  end
  integrator.qold = max(integrator.EEst,integrator.opts.qoldinit)
  return integrator.dt/q # dtnew
end

function step_accept_controller!(integrator, alg, q)
  if ispredictive(alg)
    predictive_step_accept_controller!(integrator, alg, q)
  elseif isstandard(alg)
    standard_step_accept_controller!(integrator, alg, q)
  else
    PI_step_accept_controller!(integrator, alg, q)
  end
end

function step_reject_controller!(integrator,alg)
  if ispredictive(alg)
    integrator.dt = integrator.success_iter == 0 ? 0.1integrator.dt : integrator.dt/integrator.qold
  elseif isstandard(alg)
    integrator.dt = integrator.qold
  else #PI
    integrator.dt = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
  end
  return
end

function stepsize_controller!(integrator, alg::JVODE)
  if iszero(integrator.EEst)
    η = integrator.opts.qmax
  else
    η = integrator.cache.η
    integrator.qold = η
  end
  η
end

function step_accept_controller!(integrator,alg::JVODE,η)
  q = inv(η)
  if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    q = one(q)
  end
  return integrator.dt/q  # dtnew
end

function step_reject_controller!(integrator,alg::JVODE)
  integrator.dt *= integrator.qold
end

function stepsize_controller!(integrator, alg::QNDF)
  cnt = integrator.iter
  if cnt <= 3
    q = standard_stepsize_controller!(integrator, alg)
  else
    q = integrator.dt/integrator.cache.h
    integrator.qold = integrator.dt/q
  end
  q
end


@inline function stepsize_controller!(integrator,alg::Union{ExtrapolationMidpointDeuflhard,ImplicitDeuflhardExtrapolation})
  # Dummy function
  # ExtrapolationMidpointDeuflhard's stepsize scaling is stored in the cache;
  # it is computed by  stepsize_controller_internal! (in perfom_step!) resp. stepsize_predictor!
  # (in step_accept_controller! and step_reject_controller!)
  zero(typeof(integrator.opts.qmax))
end

function stepsize_controller_internal!(integrator,alg::Union{ExtrapolationMidpointDeuflhard,ImplicitDeuflhardExtrapolation})
  # Standard stepsize controller
  # Compute and save the stepsize scaling based on the latest error estimate of the current order
  if iszero(integrator.EEst)
    q = inv(integrator.opts.qmax)
  else
    # Update gamma and beta1
    integrator.opts.beta1 = typeof(integrator.opts.beta1)(1 // (2integrator.cache.n_curr + 1))
    integrator.opts.gamma = DiffEqBase.fastpow(typeof(integrator.opts.gamma)(1 // 4),integrator.opts.beta1)
    # Compute new stepsize scaling
    qtmp = DiffEqBase.fastpow(integrator.EEst,integrator.opts.beta1) / integrator.opts.gamma
    @fastmath q = max(inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), qtmp))
  end
  integrator.cache.Q[integrator.cache.n_curr - alg.n_min + 1] = q
end

function stepsize_predictor!(integrator,alg::Union{ExtrapolationMidpointDeuflhard,ImplicitDeuflhardExtrapolation},n_new::Int)
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
    integrator.opts.gamma = DiffEqBase.fastpow(typeof(integrator.opts.gamma)(1 // 4),integrator.opts.beta1)
    # Compute new stepsize scaling
    qtmp = EEst * DiffEqBase.fastpow(DiffEqBase.fastpow(tol,(1.0 - s_curr / s_new)),integrator.opts.beta1) / integrator.opts.gamma
    @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp))
  end
  integrator.cache.Q[n_new - alg.n_min + 1] = q
end

function step_accept_controller!(integrator,alg::Union{ExtrapolationMidpointDeuflhard,ImplicitDeuflhardExtrapolation},q)
  # Compute new order and stepsize, return new stepsize
  @unpack n_min, n_max = alg
  @unpack n_curr, n_old, Q = integrator.cache
  s = integrator.cache.stage_number

  # Compute new order based on available quantities
  tmp = (n_min:n_curr) .- n_min .+ 1 # Index range of quantities computed so far
  dt_new = Vector{eltype(Q)}(undef,length(tmp)+1)
  dt_new[1:end-1] = integrator.dt ./ Q[tmp] # Store for the possible new stepsizes
  dtmin = timedepentdtmin(integrator)
  dt_new[1:end-1] = max.(dtmin, min.(abs(integrator.opts.dtmax), abs.(dt_new[1:end-1]))) # Safety scaling

  # n_new is the most efficient order of the last step
  work = s[tmp] ./ dt_new[1:end-1]
  n_new = argmin(work) + n_min - 1

  # Check if n_new may be increased
  if n_new == n_curr < min(n_max, n_old + 1) # cf. win_max in perfom_step! of the last step
    # Predict stepsize scaling for order (n_new + 1)
    stepsize_predictor!(integrator, alg, n_new+1) # Update cache.Q

    # Compute and scale the corresponding stepsize
    dt_new[end] = integrator.dt ./ Q[tmp[end]+1]
    dt_new[end] = max(dtmin, min(abs(integrator.opts.dtmax), abs.(dt_new[end])))

    # Check if (n_new  + 1) would have been more efficient than n_new
    if work[end] > s[tmp[end]+1] / dt_new[end]
      n_new = n_new + 1
    end
  end

  integrator.cache.n_curr = n_new
  dt_new[n_new - n_min + 1]
end

function step_reject_controller!(integrator, alg::Union{ExtrapolationMidpointDeuflhard,ImplicitDeuflhardExtrapolation})
  # Compute and save reduced stepsize dt_red of order n_old
  # Use the latest error estimate to predict dt_red if an estimate of order n_old is not available
  if integrator.cache.n_curr < integrator.cache.n_old
      stepsize_predictor!(integrator,alg,integrator.cache.n_old) # Update cache.Q
  end
  integrator.cache.n_curr = integrator.cache.n_old # Reset order for redoing the rejected step
  dt_red = integrator.dt / integrator.cache.Q[integrator.cache.n_old - integrator.alg.n_min + 1]
  dtmin = timedepentdtmin(integrator)
  dt_red = integrator.tdir*max(dtmin, min(abs(integrator.opts.dtmax), abs(dt_red))) # Safety scaling
  integrator.dt = dt_red
end

@inline function stepsize_controller!(integrator,alg::Union{ExtrapolationMidpointHairerWanner, ImplicitHairerWannerExtrapolation})
  # Dummy function
  # ExtrapolationMidpointHairerWanner's stepsize scaling is stored in the cache;
  # it is computed by  stepsize_controller_internal! (in perfom_step!), step_accept_controller! or step_reject_controller!
  zero(typeof(integrator.opts.qmax))
end

function stepsize_controller_internal!(integrator,alg::Union{ExtrapolationMidpointHairerWanner, ImplicitHairerWannerExtrapolation})
  # Standard stepsize controller
  # Compute and save the stepsize scaling based on the latest error estimate of the current order
  if iszero(integrator.EEst)
    q = inv(integrator.opts.qmax)
  else
    # Update gamma and beta1
    integrator.opts.beta1 = typeof(integrator.opts.beta1)(1 // (2integrator.cache.n_curr + 1))
    integrator.opts.gamma = DiffEqBase.fastpow(typeof(integrator.opts.gamma)(65 // 100),integrator.opts.beta1)
    # Compute new stepsize scaling
    qtmp = DiffEqBase.fastpow(integrator.EEst,integrator.opts.beta1) / integrator.opts.gamma
    @fastmath q = max(inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), qtmp))
  end
  integrator.cache.Q[integrator.cache.n_curr + 1] = q
end

function step_accept_controller!(integrator,alg::Union{ExtrapolationMidpointHairerWanner, ImplicitHairerWannerExtrapolation},q)
  # Compute new order and stepsize, return new stepsize
  @unpack n_min, n_max = alg
  @unpack n_curr, n_old, Q, sigma = integrator.cache
  s = integrator.cache.stage_number

  # Compute new order based on available quantities
  win_min_old = min(n_old, n_curr) - 1 # cf. win_min in perfom_step! of the last step
  tmp = win_min_old:(max(n_curr, n_old) + 1) # Index range for the new order
  dt_new = fill(zero(eltype(Q)),n_max+1)
  dt_new[tmp] = integrator.dt ./ Q[tmp] # Store for the possible new stepsizes
  dtmin = timedepentdtmin(integrator)
  dt_new[tmp] = max.(dtmin, min.(abs(integrator.opts.dtmax), abs.(dt_new[tmp]))) # Safety scaling
  work= Vector{eltype(Q)}(undef,n_max+1) # work[n] is the work for order (n-1)
  work[tmp] = s[tmp] ./ dt_new[tmp]

  # Order selection
  n_new = n_old
  if n_curr == n_min # Enforce n_min + 1 ≦ n_new
    n_new = n_min + 1
  else
    if n_curr <= n_old
      if work[n_curr-1] < sigma * work[n_curr]
        n_new = max(n_curr-1,n_old-1,n_min+1) # Enforce n_min + 1≦ n_new
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
    dt_new[n_new+1] = max(dtmin, min(abs(integrator.opts.dtmax), abs(dt_new[n_new+1])))
  end
  dt_new[n_new + 1]
end

function step_reject_controller!(integrator, alg::Union{ExtrapolationMidpointHairerWanner, ImplicitHairerWannerExtrapolation})
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
  dtmin = timedepentdtmin(integrator)
  dt_red = integrator.tdir*max(dtmin, min(abs(integrator.opts.dtmax), abs(dt_red))) # Safety scaling
  integrator.dt = dt_red
end
