@inline function initialize!(integrator,cache::LinearImplicitEulerConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@inline function perform_step!(integrator,cache::LinearImplicitEulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator

  L = update_coefficients(integrator.f,t+dt,u)

  if typeof(uprev) <: AbstractArray
    W = I - dt*L
  else
    W = 1 - dt*L
  end

  # if Bs
  k = uprev

  u = W\k

  if integrator.opts.adaptive && integrator.iter > 1
    # Use 2rd divided differences a la SPICE and Shampine
    uprev2 = integrator.uprev2
    tprev = integrator.tprev
    DD3 = ((u - uprev)/((dt)*(t+dt-tprev)) + (uprev-uprev2)/((t-tprev)*(t+dt-tprev)))
    dEst = (dt^2)*abs(DD3/6)
    integrator.EEst = dEst/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
  else
    integrator.EEst = 1
  end

  integrator.fsallast = f(t+dt,u)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end#

@inline function initialize!(integrator,cache::LinearImplicitEulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

@inline function perform_step!(integrator,cache::LinearImplicitEulerCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack W,k = cache
  mass_matrix = integrator.sol.prob.mass_matrix

  L = integrator.f
  update_coefficients!(L,t+dt,u)

  # Check is_constant before redoing
  for j in 1:length(u), i in 1:length(u)
      @inbounds W[i,j] = @muladd mass_matrix[i,j]-dt*L[i,j]
  end

  k .= uprev # + B

  integrator.alg.linsolve(vec(u),W,vec(k),true)

  if integrator.opts.adaptive && integrator.iter > 1
    # Use 2rd divided differences a la SPICE and Shampine
    uprev2 = integrator.uprev2
    tprev = integrator.tprev
    dt1 = (dt)*(t+dt-tprev)
    dt2 = (t-tprev)*(t+dt-tprev)
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds DD3 = (u[i] - uprev[i])/dt1 + (uprev[i]-uprev2[i])/dt2
      dEst = (dt^2)*abs(DD3)/6
      @inbounds k[i] = dEst/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(k)
  else
    integrator.EEst = 1
  end

  f(t+dt,u,integrator.fsallast)
  @pack integrator = t,dt,u
end
