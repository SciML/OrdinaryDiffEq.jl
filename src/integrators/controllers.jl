@inline function predictive_stepsize_controller!(integrator, alg)
  # Gustafsson predictive stepsize controller
  if iszero(integrator.EEst)
    q = inv(integrator.opts.qmax)
  else
    gamma = integrator.opts.gamma
    if typeof(alg) <: Union{RKC,IRKC,SERK2}
      fac = gamma
    else
      if alg isa Union{RadauIIA3, RadauIIA5}
        @unpack iter = integrator.cache
        @unpack maxiters = alg
      else
        @unpack iter, maxiters = integrator.cache.nlsolver
      end
      fac = min(gamma,(1+2*maxiters)*gamma/(iter+2*maxiters))
    end
    expo = 1/(get_current_adaptive_order(integrator.alg,integrator.cache) + 1)
    qtmp = DiffEqBase.fastpow(integrator.EEst,expo)/fac
    @fastmath q = DiffEqBase.value(max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp)))
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
    @fastmath q = DiffEqBase.value(max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp)))
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
    q11 = DiffEqBase.value(DiffEqBase.fastpow(EEst,beta1))
    q = q11/DiffEqBase.fastpow(qold,beta2)
    integrator.q11 = q11
    @fastmath q = DiffEqBase.value(max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),q/integrator.opts.gamma)))
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
    qacc= DiffEqBase.value(max(q,qgus))
  else
    qacc = DiffEqBase.value(q)
  end
  if integrator.opts.qsteady_min <= qacc <= integrator.opts.qsteady_max
    qacc = DiffEqBase.value(one(qacc))
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
  if iszero(integrator.EEst)
    q = inv(integrator.opts.qmax)
  else
    EEst1 = integrator.cache.EEst1
    EEst2 = integrator.cache.EEst2
    if integrator.cache.nconsteps < integrator.cache.order + 2
      integrator.cache.nconsteps += 1
      q = one(integrator.qold) #quasi-contsant steps
    else
      prev_order = integrator.cache.order
      dtnew = QNDF_stepsize_and_order!(integrator.cache, integrator.EEst, EEst1, EEst2, integrator.dt, integrator.cache.order, integrator.opts.gamma)

      if(integrator.dt != dtnew || prev_order !=integrator.cache.order)
        integrator.cache.nconsteps = 1
      end
      q = integrator.dt/dtnew
      integrator.qold = q
    end
  end
  q
end

function step_accept_controller!(integrator,alg::QNDF,q)
  #step is accepted, reset count of consecutive failed steps
  integrator.cache.consfailcnt = 0
  if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    q = one(q)
  end
  return integrator.dt/q  # dtnew
end

function step_reject_controller!(integrator,alg::QNDF)
  #append no. of consecutive failed steps
  k = integrator.cache.order
  h = integrator.dt
  integrator.cache.consfailcnt += 1
  if integrator.cache.consfailcnt > 1  # postfail
    dt_optim_failed = h/2
    integrator.cache.order = k
  end
  zₛ = 1/(5//6)  # equivalent to intergrator.opts.gamma
  zᵤ = 0.1
  Fᵤ = 10
  expo = 1/(k+1)
  z = zₛ * ((integrator.EEst)^expo)
  F = inv(z)
  if z <= 10
    hₖ = F * h
  elseif z > 10
    hₖ = 0.1 * h
  end
  hₙ = hₖ
  kₙ = k
  if k > 1
    expo = 1/k
    zₖ₋₁ = 1.3 * ((integrator.cache.EEst1)^expo)
    Fₖ₋₁ = inv(zₖ₋₁)
    if zₖ₋₁ <= 10
      hₖ₋₁ = Fₖ₋₁ * h
    elseif zₖ₋₁ > 10
      hₖ₋₁ = 0.1 * h
    end
    if hₖ₋₁ > hₖ
      hₙ = min(h,hₖ₋₁)
      kₙ = max(k-1,1)
    end
  end
  
  integrator.dt = hₙ
  integrator.cache.order = kₙ
end


# this stepsize and order controller is taken from
# Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s paper
# E. Alberdi Celaya, J. J. Anza Aguirrezabala, and P. Chatzipantelidis
function QNDF_stepsize_and_order!(cache, est, estₖ₋₁, estₖ₊₁, h, k, gamma)

  #=
  ##########################################
  dt_optim_success = h
  dt_optim_failed = h
  zₛ = 1.2
  zᵤ = 0.1
  Fᵤ = 10

  expo = 1/(k+1)
  z = zₛ * ((est)^expo)
  F = inv(z)

  hₖ₋₁ = 0.0
  hₖ₊₁ = 0.0

  # step is successful
  # precalculations
  # calculating for successful step
  if z <= zₛ=#


  
  zₛ = 1/(5//6) # equivalent to intergrator.opts.gamma
  zᵤ = 0.1
  Fᵤ = 10
  expo = 1/(k+1)
  z = zₛ * ((est)^expo)
  F = inv(z)
  hₙ = h
  if z <= zᵤ
    hₖ = Fᵤ * h
  elseif zᵤ < z <= zₛ
    hₖ = F * h
  end
  hₖ₋₁ = 0.0
  hₖ₊₁ = 0.0
  

  if k > 1
    expo = 1/k
    zₖ₋₁ = 1.3 * ((estₖ₋₁)^expo)
    Fₖ₋₁ = inv(zₖ₋₁)
    if zₖ₋₁ <= 0.1
      hₖ₋₁ = 10 * h
    elseif 0.1 < zₖ₋₁
      hₖ₋₁ = Fₖ₋₁ * h
    end
    if hₖ₋₁ > hₖ
      hₙ = hₖ₋₁
      kₙ = k-1
    end
  end

  if k < cache.max_order && cache.nconsteps > k+1
    expo = 1/(k+2)
    zₖ₊₁ = 1.4 * ((estₖ₊₁)^expo)
    Fₖ₊₁ = inv(zₖ₊₁)

    if zₖ₊₁<= 0.1
      hₖ₊₁ = 10 * h
    elseif 0.1 < zₖ₊₁
      hₖ₊₁ = Fₖ₊₁ * h
    end
    if hₖ₊₁ > hₙ
      hₙ = hₖ₊₁
      kₙ = k+1
    end
  end
  # adp order and step conditions

  if hₙ <= h
    hₙ = h
    kₙ = k
  end

  cache.order = kₙ

  return hₙ
  #=
  else
    # fail step calcuator

    # step is not successful
    if cache.consfailcnt >= 1  # postfail
      dt_optim_failed = h/2
      cache.order = k
    end

    if 1.2 < z <= 10
      hₖ = F * h
    elseif z > 10
      hₖ = 0.1 * h
    end
    hₙ = hₖ
    kₙ = k
    if k > 1
      expo = 1/k
      zₖ₋₁ = 1.3 * ((estₖ₋₁)^expo)
      Fₖ₋₁ = inv(zₖ₋₁)
      if 1.3 < zₖ₋₁ <= 10
        hₖ₋₁ = Fₖ₋₁ * h
      elseif zₖ₋₁ > 10
        hₖ₋₁ = 0.1 * h
      end

      if hₖ₋₁ > hₖ
        hₙ = min(h,hₖ₋₁)
        kₙ = max(k-1,1)
      end
    end
    dt_optim_failed = hₙ
    cache.order = kₙ
  end
  return dt_optim_success, dt_optim_failed
  =#
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

@inline function stepsize_controller!(integrator,alg::Union{ExtrapolationMidpointHairerWanner, ImplicitHairerWannerExtrapolation, ImplicitEulerExtrapolation, ImplicitEulerBarycentricExtrapolation})
  # Dummy function
  # ExtrapolationMidpointHairerWanner's stepsize scaling is stored in the cache;
  # it is computed by  stepsize_controller_internal! (in perfom_step!), step_accept_controller! or step_reject_controller!
  zero(typeof(integrator.opts.qmax))
end

function stepsize_controller_internal!(integrator,alg::Union{ExtrapolationMidpointHairerWanner, ImplicitHairerWannerExtrapolation, ImplicitEulerExtrapolation, ImplicitEulerBarycentricExtrapolation})
  # Standard stepsize controller
  # Compute and save the stepsize scaling based on the latest error estimate of the current order
  if typeof(alg) <: Union{ImplicitEulerExtrapolation,ImplicitEulerBarycentricExtrapolation,ImplicitHairerWannerExtrapolation}
    if iszero(integrator.EEst)
      q = inv(integrator.opts.qmax)
    else
      # Update gamma and beta1
      if typeof(alg) <: ImplicitHairerWannerExtrapolation
        integrator.opts.beta1 = typeof(integrator.opts.beta1)(1 // (2integrator.cache.n_curr + 1))
      elseif typeof(alg) <: ImplicitEulerExtrapolation
        integrator.opts.beta1 = typeof(integrator.opts.beta1)(1 // (integrator.cache.n_curr))
      else
        integrator.opts.beta1 = typeof(integrator.opts.beta1)(1 // (integrator.cache.n_curr - 1))
      end
      integrator.opts.gamma = DiffEqBase.fastpow(typeof(integrator.opts.gamma)(65 // 100),integrator.opts.beta1)
      # Compute new stepsize scaling
      qtmp = DiffEqBase.fastpow(integrator.EEst,integrator.opts.beta1) / (integrator.opts.gamma)
      @fastmath q = max(inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), qtmp))
    end
    integrator.cache.Q[integrator.cache.n_curr + 1] = q
  else
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
end

function step_accept_controller!(integrator,alg::Union{ExtrapolationMidpointHairerWanner, ImplicitHairerWannerExtrapolation, ImplicitEulerExtrapolation, ImplicitEulerBarycentricExtrapolation},q)
  # Compute new order and stepsize, return new stepsize
  @unpack n_min, n_max = alg
  @unpack n_curr, n_old, Q, sigma, work, dt_new = integrator.cache
  s = integrator.cache.stage_number

  # Compute new order based on available quantities
  win_min_old = min(n_old, n_curr) - 1 # cf. win_min in perfom_step! of the last step
  tmp = win_min_old:(max(n_curr, n_old) + 1) # Index range for the new order
  #@show size(dt_new)
  fill!(dt_new, zero(eltype(dt_new)))
  @.. Q = integrator.dt/Q
  copyto!(dt_new,win_min_old,Q,win_min_old,(max(n_curr, n_old) + 1) - win_min_old + 1)
  @.. Q = integrator.dt/Q
  dtmin = timedepentdtmin(integrator)
  fill!(work,zero(eltype(work))) # work[n] is the work for order (n-1)
  for i in tmp
    work[i] = s[i]/dt_new[i]
  end
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

function step_reject_controller!(integrator, alg::Union{ExtrapolationMidpointHairerWanner, ImplicitHairerWannerExtrapolation, ImplicitEulerExtrapolation, ImplicitEulerBarycentricExtrapolation})
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
