function initialize!(integrator,cache::RichardsonEulerCache)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # For the interpolation, needs k at the updated point

  cache.step_no = 1
  cache.cur_order = max(integrator.alg.init_order, integrator.alg.min_order)
end

function perform_step!(integrator,cache::RichardsonEulerCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,T,u_tilde,atmp,dtpropose,cur_order,A = cache

  @muladd @. u = uprev + dt*fsalfirst

  T[1,1] = copy(u)

  halfdt = dt/2
  for i in 2:size(T)[1]
    # Solve using Euler method
    @muladd @. u = uprev + halfdt*fsalfirst
    f(k, u, p, t+halfdt)
    for j in 2:2^(i-1)
      @muladd @. u = u + halfdt*k
      f(k, u, p, t+j*halfdt)
    end
    T[i,1] = copy(u)
    # Richardson Extrapolation
    for j in 2:i
      T[i,j] = ((2^(j-1))*T[i,j-1] - T[i-1,j-1])/((2^(j-1)) - 1)
    end

    halfdt = halfdt/2
  end

  if integrator.opts.adaptive
      minimum_work = Inf
      range_start = max(2,cur_order - 1)
      if cache.step_no == one(cache.step_no)
          range_start = 2
      end

      for i = range_start:min(size(T)[1], cur_order + 1)

          A = 2^(i-1)
          @. u_tilde = T[i,i] - T[i,i-1]
          atmp = calculate_residuals(u_tilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
          EEst = integrator.opts.internalnorm(atmp)

          beta1 = integrator.opts.beta1
          e = integrator.EEst
          qold = integrator.qold

          integrator.opts.beta1 = 1/(i+1)
          integrator.EEst = EEst
          dtpropose = step_accept_controller!(integrator,integrator.alg,stepsize_controller!(integrator,integrator.alg))
          integrator.EEst = e
          integrator.opts.beta1 = beta1
          integrator.qold = qold

          work = A/dtpropose

          if work < minimum_work
              integrator.opts.beta1 = 1/(i+1)
              cache.dtpropose = dtpropose
              cache.cur_order = i
              minimum_work = work
              integrator.EEst = EEst
          end
      end
  end

  # using extrapolated value of u
  @. u = T[cache.cur_order, cache.cur_order]
  cache.step_no = cache.step_no + 1
  f(k, u, p, t+dt)
end

function initialize!(integrator,cache::RichardsonEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  cache.step_no = 1
  cache.cur_order = max(integrator.alg.init_order, integrator.alg.min_order)
end

function perform_step!(integrator,cache::RichardsonEulerConstantCache,repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  @unpack dtpropose, T, cur_order, work, A = cache
  @muladd u = @. uprev + dt*integrator.fsalfirst
  T[1,1] = u
  halfdt = dt/2
  for i in 2:min(size(T)[1], cur_order+1)
    # Solve using Euler method
    @muladd u = @. uprev + halfdt*integrator.fsalfirst
    k = f(u, p, t+halfdt)

    for j in 2:2^(i-1)
      @muladd u = @. u + halfdt*k
      k = f(u, p, t+j*halfdt)
    end
    T[i,1] = u

    # Richardson Extrapolation
    for j in 2:i
      T[i,j] = ((2^(j-1))*T[i,j-1] - T[i-1,j-1])/((2^(j-1)) - 1)
    end
    halfdt = halfdt/2
  end

  if integrator.opts.adaptive
      minimum_work = Inf
      range_start = max(2,cur_order - 1)
      if cache.step_no == one(cache.step_no)
          range_start = 2
      end

      for i = range_start:min(size(T)[1], cur_order + 1)

          A = 2^(i-1)
          u_tilde = T[i,i] - T[i,i-1]
          atmp = calculate_residuals(u_tilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
          EEst = integrator.opts.internalnorm(atmp)

          beta1 = integrator.opts.beta1
          e = integrator.EEst
          qold = integrator.qold

          integrator.opts.beta1 = 1/(i+1)
          integrator.EEst = EEst
          dtpropose = step_accept_controller!(integrator,integrator.alg,stepsize_controller!(integrator,integrator.alg))
          integrator.EEst = e
          integrator.opts.beta1 = beta1
          integrator.qold = qold

          work = A/dtpropose

          if work < minimum_work
              integrator.opts.beta1 = 1/(i+1)
              cache.dtpropose = dtpropose
              cache.cur_order = i
              minimum_work = work
              integrator.EEst = EEst
          end
      end
  end

  cache.step_no = cache.step_no + 1
  integrator.u = T[cache.cur_order,cache.cur_order]

  k = f(integrator.u, p, t+dt)
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast

  # Use extrapolated value of u
end

function initialize!(integrator,cache::ExtrapolationMidpointDeuflhardCache)
  # cf. initialize! of MidpointCache
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

function perform_step!(integrator, cache::ExtrapolationMidpointDeuflhardCache, repeat_step = false)
  @unpack t, uprev, dt, f, p = integrator
  @unpack n_curr  = cache
  tol = integrator.opts.internalnorm(integrator.opts.reltol, t) # Deuflhard's approach relies on EEstD ≈ ||relTol||
  fill!(cache.T,zero(uprev))


  if integrator.opts.adaptive
    # Set up the order window
    win_min = max(integrator.alg.n_min, n_curr - 1)
    win_max = min(integrator.alg.n_max, n_curr + 1)

    # Set up the current extrapolation order
    cache.n_old = n_curr # Save the suggested order for step_accept_controller! and step_reject_controller!
    n_curr = win_min # Start with smallest order in the order window
  end

  #Compute the internal discretisations
  for i = 0 : n_curr
    j_int = 2Int64(cache.subdividing_sequence[i+1])
    dt_int = dt / (2j_int) # stepsize of the ith internal discretisation
     cache.u_temp2 = uprev
     cache.u_temp1 = cache.u_temp2 + dt_int*cache.fsalfirst # Euler starting step
    for j = 2 : 2j_int
    f(cache.k, cache.u_temp1, p, t + (j-1)dt_int)
      cache.T[i+1] = cache.u_temp2 + 2dt_int*cache.k # Explicit Midpoint rule
       cache.u_temp2 = cache.u_temp1
       cache.u_temp1 = cache.T[i+1]
    end
  end

  if integrator.opts.adaptive
    # Compute all information relating to an extrapolation order ≦ win_min
    for i = integrator.alg.n_min:n_curr
      integrator.u = eltype(uprev).(cache.extrapolation_scalars[i+1]) * sum( broadcast(*, cache.T[1:(i+1)], eltype(uprev).(cache.extrapolation_weights[1:(i+1), (i+1)])) ) # Approximation of extrapolation order i
      cache.u_tilde = eltype(uprev).(cache.extrapolation_scalars_2[i]) * sum( broadcast(*, cache.T[2:(i+1)], eltype(uprev).(cache.extrapolation_weights_2[1:i, i])) ) # and its internal counterpart
      calculate_residuals!(cache.res, integrator.u, cache.u_tilde, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      integrator.EEst = integrator.opts.internalnorm(cache.res, t)
      cache.n_curr = i # Update chache's n_curr for stepsize_controller_internal!
      stepsize_controller_internal!(integrator, integrator.alg) # Update cache.Q
    end

    # Check if an approximation of some order in the order window can be accepted
    while n_curr <= win_max
      if integrator.EEst <= 1.0
        # Accept current approximation u of order n_curr
        break
      elseif integrator.EEst <= tol^(cache.stage_number[n_curr - integrator.alg.n_min + 1] / cache.stage_number[win_max - integrator.alg.n_min + 1] - 1)
        # Reject current approximation order but pass convergence monitor
        # Compute approximation of order (n_curr + 1)
        n_curr = n_curr + 1
        cache.n_curr = n_curr

        # Update cache.T
        j_int = 2Int64(cache.subdividing_sequence[n_curr + 1])
        dt_int = dt / (2j_int) # stepsize of the new internal discretisation
         cache.u_temp2 = uprev
         cache.u_temp1 = cache.u_temp2 + dt_int * cache.fsalfirst # Euler starting step
        for j = 2 : 2j_int
          f(cache.k, cache.u_temp1, p, t + (j-1)dt_int)
          cache.T[n_curr+1] = cache.u_temp2 + 2dt_int*cache.k
           cache.u_temp2 = cache.u_temp1
           cache.u_temp1 = cache.T[n_curr+1]
        end

        # Update u, integrator.EEst and cache.Q
        integrator.u = eltype(uprev).(cache.extrapolation_scalars[n_curr+1]) * sum( broadcast(*, cache.T[1:(n_curr+1)], eltype(uprev).(cache.extrapolation_weights[1:(n_curr+1), (n_curr+1)])) ) # Approximation of extrapolation order n_curr
        cache.u_tilde = eltype(uprev).(cache.extrapolation_scalars_2[n_curr]) * sum( broadcast(*, cache.T[2:(n_curr+1)], eltype(uprev).(cache.extrapolation_weights_2[1:n_curr, n_curr])) ) # and its internal counterpart
        calculate_residuals!(cache.res, integrator.u, cache.u_tilde, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(cache.res, t)
        stepsize_controller_internal!(integrator, integrator.alg) # Update cache.Q
      else
          # Reject the current approximation and not pass convergence monitor
          break
      end
    end
  else
    integrator.u = eltype(uprev).(cache.extrapolation_scalars[n_curr+1]) * sum( broadcast(*, cache.T[1:(n_curr+1)], eltype(uprev).(cache.extrapolation_weights[1:(n_curr+1), (n_curr+1)])) ) # Approximation of extrapolation order n_curr
  end

  f(cache.k, integrator.u, p, t+dt) # Update FSAL
end

function initialize!(integrator,cache::ExtrapolationMidpointDeuflhardConstantCache)
  # cf. initialize! of MidpointConstantCache
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::ExtrapolationMidpointDeuflhardConstantCache, repeat_step=false)
  @unpack t, uprev, dt, f, p = integrator
  @unpack n_curr = cache
  T = fill(zero(uprev), integrator.alg.n_max + 1) # Storage for the internal discretisations obtained by the explicit midpoint rule
  u_temp1, u_temp2 = copy(uprev), copy(uprev) # Auxiliary variables for computing the internal discretisations
  u, u_tilde = copy(uprev), copy(uprev) # Storage for the latest approximation and its internal counterpart
  tol = integrator.opts.internalnorm(integrator.opts.reltol, t) # Deuflhard's approach relies on EEstD ≈ ||relTol||


  if integrator.opts.adaptive
    # Set up the order window
    win_min = max(integrator.alg.n_min, n_curr - 1)
    win_max = min(integrator.alg.n_max, n_curr + 1)

    # Set up the current extrapolation order
    cache.n_old = n_curr # Save the suggested order for step_accept_controller! and step_reject_controller!
    n_curr = win_min # Start with smallest order in the order window
  end

  # Compute the internal discretisations
  for i = 0 : n_curr
    j_int = 2Int64(cache.subdividing_sequence[i+1])
    dt_int = dt / (2j_int) # stepsize of the ith internal discretisation
    u_temp2 = uprev
    u_temp1 = u_temp2 + dt_int*integrator.fsalfirst # Euler starting step
    for j = 2 : 2j_int
      T[i+1] = u_temp2 + 2dt_int*f(u_temp1, p, t + (j-1)dt_int) # Explicit Midpoint rule
      u_temp2 = u_temp1
      u_temp1 = T[i+1]
    end
  end

  if integrator.opts.adaptive
    # Compute all information relating to an extrapolation order ≦ win_min
    for i = integrator.alg.n_min:n_curr
      u = eltype(uprev).(cache.extrapolation_scalars[i+1]) * sum( broadcast(*, T[1:(i+1)], eltype(uprev).(cache.extrapolation_weights[1:(i+1), (i+1)])) ) # Approximation of extrapolation order i
      u_tilde = eltype(uprev).(cache.extrapolation_scalars_2[i]) * sum( broadcast(*, T[2:(i+1)], eltype(uprev).(cache.extrapolation_weights_2[1:i, i])) ) # and its internal counterpart
      res = calculate_residuals(u, u_tilde, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
      integrator.EEst = integrator.opts.internalnorm(res, t)
      cache.n_curr = i # Update chache's n_curr for stepsize_controller_internal!
      stepsize_controller_internal!(integrator, integrator.alg) # Update cache.Q
    end

    # Check if an approximation of some order in the order window can be accepted
    while n_curr <= win_max
      if integrator.EEst <= 1.0
        # Accept current approximation u of order n_curr
        break
      elseif integrator.EEst <= tol^(cache.stage_number[n_curr - integrator.alg.n_min + 1] / cache.stage_number[win_max - integrator.alg.n_min + 1] - 1)
        # Reject current approximation order but pass convergence monitor
        # Compute approximation of order (n_curr + 1)
        n_curr = n_curr + 1
        cache.n_curr = n_curr

        # Update T
        j_int = 2Int64(cache.subdividing_sequence[n_curr + 1])
        dt_int = dt / (2j_int) # stepsize of the new internal discretisation
        u_temp2 = uprev
        u_temp1 = u_temp2 + dt_int * integrator.fsalfirst # Euler starting step
        for j = 2 : 2j_int
          T[n_curr+1] = u_temp2 + 2dt_int*f(u_temp1, p, t + (j-1)dt_int)
          u_temp2 = u_temp1
          u_temp1 = T[n_curr+1]
        end

        # Update u, integrator.EEst and cache.Q
        u = eltype(uprev).(cache.extrapolation_scalars[n_curr+1]) * sum( broadcast(*, T[1:(n_curr+1)], eltype(uprev).(cache.extrapolation_weights[1:(n_curr+1), (n_curr+1)])) ) # Approximation of extrapolation order n_curr
        u_tilde = eltype(uprev).(cache.extrapolation_scalars_2[n_curr]) * sum( broadcast(*, T[2:(n_curr+1)], eltype(uprev).(cache.extrapolation_weights_2[1:n_curr, n_curr])) ) # and its internal counterpart
        res = calculate_residuals(u, u_tilde, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(res, t)
        stepsize_controller_internal!(integrator, integrator.alg) # Update cache.Q
      else
          # Reject the current approximation and not pass convergence monitor
          break
      end
    end
  else
    u = eltype(uprev).(cache.extrapolation_scalars[n_curr+1]) * sum( broadcast(*, T[1:(n_curr+1)], eltype(uprev).(cache.extrapolation_weights[1:(n_curr+1), (n_curr+1)])) ) # Approximation of extrapolation order n_curr
  end


  # Save the latest approximation and update FSAL
  integrator.u = u
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end
