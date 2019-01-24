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

  @. u = uprev + dt*fsalfirst
  @. T[1,1] = u

  halfdt = dt/2
  for i in 2:m
    # Solve using Euler method
    @. u = uprev + halfdt*fsalfirst
    f(k, u, p, t+halfdt)

    for j in 2:2^(i-1)
      @. u = u + halfdt*k
      f(k, u, p, t+j*halfdt)
    end
    @. T[i,1] = u

    # Richardson Extrapolation
    for j in 2:i
      @. T[i,j] = (2^(2*i)*T[i,j-1] - T[i-1,j-1])/(2^(2*i) - 1)
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
end

function perform_step!(integrator,cache::RichardsonEulerConstantCache,repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  @unpack m = cache
  @muladd u = @. uprev + dt*integrator.fsalfirst
  T = zeros(typeof(uprev), (m,m))
  T[1,1] = u

  halfdt = dt/2
  for i in 2:m
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
      T[i,j] = (2^(2*i)*T[i,j-1] - T[i-1,j-1])/(2^(2*i) - 1)
    end
    halfdt = halfdt/2
  end

  k = f(T[m,m], p, t+dt)

  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  # Use extrapolated value of u
  integrator.u = T[m,m]
end
