function initialize!(integrator, cache::MagnusMidpointCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

function perform_step!(integrator, cache::MagnusMidpointCache, repeat_step=false)
  @unpack t,dt,uprev,u,p,alg = integrator
  @unpack W,k,tmp = cache
  mass_matrix = integrator.f.mass_matrix

  L = integrator.f.f
  update_coefficients!(L,u,p,t+dt/2)

  if integrator.alg.krylov
    u .= expv(dt, L, u; m=min(alg.m, size(L,1)), opnorm=integrator.opts.internalopnorm, iop=alg.iop)
  else
    A = Matrix(L) #size(L) == () ? convert(Number, L) : convert(AbstractMatrix, L)
    u .= exp(dt*L) * u
  end

  integrator.f(integrator.fsallast,u,p,t+dt)
  integrator.destats.nf += 1
end

function initialize!(integrator, cache::MagnusLeapfrogCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

function perform_step!(integrator, cache::MagnusLeapfrogCache, repeat_step=false, alg_extrapolates=true, iter=1)
  @unpack t,dt,uprev,uprev2,u,p,alg,iter = integrator
  @unpack W,k,tmp = cache
  mass_matrix = integrator.f.mass_matrix
    # println("iter   : $iter")
  if iter==1
    L = integrator.f.f
    update_coefficients!(L,u,p,t+dt/2)
    if integrator.alg.krylov
      u .= expv(dt, L, u; m=min(alg.m, size(L,1)), opnorm=integrator.opts.internalopnorm, iop=alg.iop)
    else
      A = Matrix(L) #size(L) == () ? convert(Number, L) : convert(AbstractMatrix, L)
      u .= exp(dt*L) * u
    end
  
    integrator.f(integrator.fsallast,u,p,t+dt)
    integrator.destats.nf += 1
    iter += 1
  else
    L = integrator.f.f
    update_coefficients!(L,u,p,t)
    if integrator.alg.krylov
      u .= expv(2*dt, L, uprev2; m=min(alg.m, size(L,1)), opnorm=integrator.opts.internalopnorm, iop=alg.iop)
    else
      A = Matrix(L) #size(L) == () ? convert(Number, L) : convert(AbstractMatrix, L)
      u .= exp(2*dt*L) * uprev2
    end
    uprev=u
    integrator.f(integrator.fsallast,u,p,t+dt)
    integrator.destats.nf += 1
  end
end

function initialize!(integrator, cache::LinearExponentialConstantCache)
  # Pre-start fsal
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
  integrator.fsallast = zero(integrator.fsalfirst)

  # Initialize interpolation derivatives
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::LinearExponentialConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  alg = unwrap_alg(integrator, true)
  A = f.f # assume f to be an ODEFunction wrapped around a linear operator

  if alg.krylov == :off
    u = exp(dt * Matrix(f)) * integrator.u
  elseif alg.krylov == :simple
    u = expv(dt, A, integrator.u; m=min(alg.m, size(A,1)), opnorm=integrator.opts.internalopnorm, iop=alg.iop)
  else
    u = expv_timestep(dt, A, integrator.u; m=min(alg.m, size(A,1)), iop=alg.iop,
                      opnorm=integrator.opts.internalopnorm, tol=integrator.opts.reltol)
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::LinearExponentialCache)
  # Pre-start fsal
  integrator.fsalfirst = zero(cache.rtmp)
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
  integrator.fsallast = zero(integrator.fsalfirst)

  # Initialize interpolation derivatives
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::LinearExponentialCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp, KsCache = cache
  alg = unwrap_alg(integrator, true)
  A = f.f # assume f to be an ODEFunction wrapped around a linear operator

  if alg.krylov == :off
    E = exp(dt * Matrix(A))
    mul!(tmp, E, u)
  elseif alg.krylov == :simple
    Ks, expv_cache = KsCache
    arnoldi!(Ks, A, u; m=min(alg.m, size(A,1)), opnorm=integrator.opts.internalopnorm, iop=alg.iop)
    expv!(tmp, dt, Ks; cache=expv_cache)
  else
    expv_timestep!(tmp, dt, A, u; adaptive=true, caches=KsCache, m=min(alg.m, size(A,1)), iop=alg.iop,
                   opnorm=integrator.opts.internalopnorm, tol=integrator.opts.reltol)
  end

  # Update integrator state
  u .= tmp
  f(integrator.fsallast, u, p, t + dt)
  integrator.destats.nf += 1
  # integrator.k is automatically set due to aliasing
end
