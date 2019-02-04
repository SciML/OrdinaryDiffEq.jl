function initialize!(integrator,cache::RichardsonEulerCache)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # For the interpolation, needs k at the updated point
end

function perform_step!(integrator,cache::RichardsonEulerCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,T = cache
  @unpack m = cache.tab

  @muladd @. u = uprev + dt*fsalfirst
  T[1,1] = copy(u)
  halfdt = dt/2
  for i in 2:m
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

  # using extrapolated value of u
  @. u = T[m,m]

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

  cache.cur_order = integrator.alg.init_order
end

function perform_step!(integrator,cache::RichardsonEulerConstantCache,repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  @unpack dtpropose, T, cur_order, prev_work, work, A = cache
  @muladd u = @. uprev + dt*integrator.fsalfirst
  A = 1 + 1

  T[1,1] = u
  dtpropose = dt #doubt but error_1 isn't defined
  work = A/dtpropose

  halfdt = dt/2
  for i in 2:size(T)[1]
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

    integrator.EEst = abs(T[i,i] - T[i,i-1])
    prevdtpropose = dtpropose
    dtpropose = dt*0.94*(0.65/integrator.EEst)^(1/(2*i - 1))
    A = A + dt/halfdt
    prev_work = work
    work = A/dtpropose

    if integrator.opts.adaptive && cur_order > 2
      if i == cur_order - 1
          if integrator.EEst <= one(integrator.EEst)
              if work < 0.9*prev_work
                  cache.dtpropose = dtpropose*((A + (2*dt/halfdt))/A)
              else
                  cache.cur_order = cur_order - 1
                  cache.dtpropose = dtpropose
              end
              integrator.u = T[i,i]
              k = f(T[i,i], p, t+dt)

              integrator.fsallast = k
              integrator.k[1] = integrator.fsalfirst
              integrator.k[2] = integrator.fsallast
              return
          else
              if integrator.EEst > 4*(dt/halfdt)^4
                  cache.cur_order = cur_order - 1
                  cache.dtpropose = dtpropose
                  return
              end
          end
      elseif i == cur_order
          if integrator.EEst <= one(integrator.EEst)
              if prev_work < 0.9*work
                  cache.cur_order = cur_order - 1
                  cache.dtpropose = prevdtpropose
              elseif work < 0.9*prev_work
                  cache.cur_order = min(integrator.alg.max_order, cur_order + 1)
                  cache.dtpropose = dtpropose*((A + 2*dt/halfdt)/A)
              else
                  cache.dtpropose = dtpropose
              end
              integrator.u = T[i,i]
              k = f(T[i,i], p, t+dt)

              integrator.fsallast = k
              integrator.k[1] = integrator.fsalfirst
              integrator.k[2] = integrator.fsallast
              return
          else
              if integrator.EEst > (2*dt/halfdt)^2
                  cache.dtpropose = dtpropose
                  return
              else
                  #setting variables for next condition
                  if prev_work < 0.9*work
                      work = prev_work
                      cur_order = cur_order - 1
                  end
              end
          end
      else
          if integrator.EEst <= one(integrator.EEst)
            cache.cur_order = cur_order
            if work < 0.9*prev_work
                cache.cur_order = min(integrator.alg.max_order,cache.cur_order + 1)
            end
            k = f(T[i,i], p, t+dt)

            integrator.fsallast = k
            integrator.k[1] = integrator.fsalfirst
            integrator.k[2] = integrator.fsallast

            integrator.u = T[i,i]
            return
          else
              cache.dtpropose = dtpropose
          end
      end
    end

    halfdt = halfdt/2
  end

  k = f(T[integrator.alg.max_order,integrator.alg.max_order], p, t+dt)

  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  # Use extrapolated value of u
  integrator.u = T[integrator.alg.max_order,integrator.alg.max_order]
  cache.dtpropose = dt
end
