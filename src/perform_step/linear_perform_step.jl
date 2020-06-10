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
function initialize!(integrator, cache::MagnusGauss4Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

function perform_step!(integrator, cache::MagnusGauss4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,p,alg = integrator
  @unpack W,k,tmp = cache
  mass_matrix = integrator.f.mass_matrix
  # println("**************")
  # println("Initial uprev : $uprev and Initial u : $u")
  L1 = deepcopy(integrator.f.f)
  L2 = deepcopy(integrator.f.f)
  update_coefficients!(L1,uprev,p,t+dt*(1/2+sqrt(3)/6))
  A = Matrix(L1)
  update_coefficients!(L2,uprev,p,t+dt*(1/2-sqrt(3)/6))
  B = Matrix(L2)
  if integrator.alg.krylov
    u .= expv(dt,(A+B) ./ 2 + (dt*sqrt(3)) .* (B*A-A*B) ./ 12, u; m=min(alg.m, size(L1,1)), opnorm=integrator.opts.internalopnorm, iop=alg.iop)
  else
    #A = Matrix(L) #size(L) == () ? convert(Number, L) : convert(AbstractMatrix, L)
    u .= exp((dt/2) .* (A+B)+((dt^2)*(sqrt(3)/12)) .* (B*A-A*B)) * uprev
  end
  integrator.f(integrator.fsallast,u,p,t+dt)
  integrator.destats.nf += 1
  # println("Final uprev : $uprev and Final u : $u")
  # println("**************")
end

function initialize!(integrator, cache::LieEulerCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

function perform_step!(integrator, cache::LieEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,p,alg = integrator
  @unpack W,k,tmp = cache
  mass_matrix = integrator.f.mass_matrix

  L = integrator.f.f
  update_coefficients!(L,u,p,t)

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

cay!(tmp, A) = mul!(tmp, inv(I - 1/2 * A), (I + 1/2 * A))
cay(A) = inv(I - 1/2 * A) * (I + 1/2 * A)

function initialize!(integrator, cache::CayleyEulerConstantCache)
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

function perform_step!(integrator, cache::CayleyEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  alg = unwrap_alg(integrator, true)

  if f isa SplitFunction
    A = f.f1.f
  else  # f isa ODEFunction
    A = f.f
  end

  L = update_coefficients(A, uprev, p, t)
  V = cay(L*dt)
  u = V * uprev * transpose(V)

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::CayleyEulerCache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
  integrator.destats.nf += 1
end

function perform_step!(integrator, cache::CayleyEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,p,f,alg = integrator
  @unpack k,V,tmp = cache
  mass_matrix = integrator.f.mass_matrix

  if f isa SplitFunction
    L = f.f1.f
  else  # f isa ODEFunction
    L = f.f
  end

  update_coefficients!(L, uprev, p, t)

  cay!(V, L*dt)
  mul!(tmp, uprev, transpose(V))
  mul!(u, V, tmp)

  # Update integrator state
  integrator.f(integrator.fsallast,u,p,t+dt)
  integrator.destats.nf += 1
end
