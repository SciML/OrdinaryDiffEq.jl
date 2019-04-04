@inline function predictive_stepsize_controller!(integrator, alg)
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
        @unpack nl_iters, max_iter = integrator.cache.nlsolver
      end
      fac = min(gamma,(1+2*max_iter)*gamma/(nl_iters+2*max_iter))
    end
    expo = 1/(get_current_adaptive_order(integrator.alg,integrator.cache)+1)
    qtmp = (integrator.EEst^expo)/fac
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
    qtmp = integrator.EEst^(1/(get_current_adaptive_order(integrator.alg,integrator.cache)+1))/integrator.opts.gamma
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
    @fastmath q11 = EEst^beta1
    @fastmath q = q11/(qold^beta2)
    integrator.q11 = q11
    @fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),q/integrator.opts.gamma))
  end
  q
end

function stepsize_controller!(integrator,alg)
  if ispredictive(alg)
    predictive_stepsize_controller!(integrator, alg)
  elseif isstandard(alg)
    standard_stepsize_controller!(integrator, alg)
  else # Default is PI-controller
    PI_stepsize_controller!(integrator, alg)
  end
end

@inline function predictive_step_accept_controller!(integrator, alg, q)
  if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    qacc = one(q)
  end
  if integrator.success_iter > 0
    expo = 1/(get_current_adaptive_order(integrator.alg,integrator.cache)+1)
    qgus = (integrator.dtacc/integrator.dt)*(((integrator.EEst^2)/integrator.erracc)^expo)
    qgus = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qgus/integrator.opts.gamma))
    qacc = max(q,qgus)
  else
    qacc = q
  end
  integrator.dtacc = integrator.dt
  integrator.erracc = max(1e-2,integrator.EEst)
  return integrator.dt/qacc
end

standard_step_accept_controller!(integrator, alg, q) = integrator.dt/q

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
  return η * integrator.dt  # dtnew
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
