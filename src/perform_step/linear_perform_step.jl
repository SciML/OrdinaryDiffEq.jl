function initialize!(integrator, cache::MidpointSplittingCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

function perform_step!(integrator, cache::MidpointSplittingCache, repeat_step=false)
  @unpack t,dt,uprev,u = integrator
  @unpack W,k,tmp = cache
  mass_matrix = integrator.f.mass_matrix

  L = integrator.f
  update_coefficients!(L,u,p,t+dt/2)

  A = L.As[1]
  Bs = L.As[2:end]

  copyto!(tmp, uprev)
  for B in reverse(Bs)
    u .= exp((dt/2)*B)*tmp
    @swap!(tmp,u)
  end

  u .= exp(dt*A)*tmp

  for B in Bs
    tmp .= exp((dt/2)*B)*u
    @swap!(u,tmp)
  end

  f(integrator.fsallast,u,p,t+dt)
  integrator.destats.nf += 1
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