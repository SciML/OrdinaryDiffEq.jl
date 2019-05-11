# This function calculates the largest eigenvalue
# (absolute value wise) by power iteration.
const RKCAlgs = Union{RKC,IRKC,ESERK5}
function maxeig!(integrator, cache::OrdinaryDiffEqConstantCache)
  isfirst = integrator.iter == 1 || integrator.u_modified
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  maxiter = (integrator.alg isa ESERK5) ? 100 : 50

  safe = (typeof(integrator.alg) <: RKCAlgs) ? 1.0 : 1.2
  # Initial guess for eigenvector `z`
  if isfirst
    if typeof(integrator.alg) <: RKCAlgs
      if integrator.alg isa IRKC
        z = cache.du₂
      else
        z = fsalfirst
      end
    else
      fz = fsalfirst
      z = f(fz, p, t)
      integrator.destats.nf += 1
    end
  else
    z = cache.zprev
  end
  # Perturbation
  u_norm = integrator.opts.internalnorm(uprev,t)
  z_norm = integrator.opts.internalnorm(z,t)
  pert   = eps(u_norm)
  sqrt_pert = sqrt(pert)
  is_u_zero = u_norm == zero(u_norm)
  is_z_zero = z_norm == zero(z_norm)
  # Normalize `z` such that z-u lie in a circle
  if ( !is_u_zero && !is_z_zero )
    dz_u = u_norm * sqrt_pert
    quot = dz_u/z_norm
    z = uprev + quot*z
  elseif !is_u_zero
    dz_u = u_norm * sqrt_pert
    z = uprev + uprev*dz_u
  elseif !is_z_zero
    dz_u = pert
    quot = dz_u/z_norm
    z *= quot
  else
    dz_u = pert
    z = dz_u
  end # endif
  # Start power iteration
  integrator.eigen_est = 0
  for iter in 1:maxiter
    if integrator.alg isa IRKC
      fz = f.f2(z, p, t)
      tmp = fz - cache.du₂
    else
      fz = f(z, p, t)
      integrator.destats.nf += 1
      tmp = fz - fsalfirst
    end
    Δ  = integrator.opts.internalnorm(tmp,t)
    eig_prev = integrator.eigen_est
    integrator.eigen_est = Δ/dz_u * safe
    # Convergence
    if typeof(integrator.alg) <: RKCAlgs # To match the constants given in the paper
      if iter >= 2 && abs(eig_prev - integrator.eigen_est) < max(integrator.eigen_est,1.0/integrator.opts.dtmax)*0.01
        integrator.eigen_est *= 1.2
        # Store the eigenvector
        cache.zprev = z - uprev
        return true
      end
    else
      if iter >= 2 && abs(eig_prev - integrator.eigen_est) < integrator.eigen_est*0.05
        # Store the eigenvector
        cache.zprev = z
        return true
      end
    end

    # Next `z`
    if Δ != zero(Δ)
      quot = dz_u/Δ
      z = uprev + quot*tmp
    else
      # An arbitrary change on `z`
      nind = length(z)
      if (nind != 1)
        ind = 1 + iter % nind
        z[ind] = uprev[ind] - (z[ind] - uprev[ind])
      else
        z = -z
      end
    end
  end
  return false
end

function maxeig!(integrator, cache::OrdinaryDiffEqMutableCache)
  isfirst = integrator.iter == 1 || integrator.u_modified
  @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
  fz, z, atmp = cache.k, cache.tmp, cache.atmp
  ccache = cache.constantcache
  maxiter = (integrator.alg isa ESERK5) ? 100 : 50
  safe = (typeof(integrator.alg) <: RKCAlgs) ? 1.0 : 1.2
  # Initial guess for eigenvector `z`
  if isfirst
    if typeof(integrator.alg) <: RKCAlgs
      if integrator.alg isa IRKC
        @.. z = cache.du₂
      else
        @.. z = fsalfirst
      end
    else
      @.. fz = u
      f(z, fz, p, t)
      integrator.destats.nf += 1
    end
  else
    @.. z = ccache.zprev
  end
  # Perturbation
  u_norm = integrator.opts.internalnorm(uprev,t)
  z_norm = integrator.opts.internalnorm(z,t)
  pert   = eps(u_norm)
  sqrt_pert = sqrt(pert)
  is_u_zero = u_norm == zero(u_norm)
  is_z_zero = z_norm == zero(z_norm)
  # Normalize `z` such that z-u lie in a circle
  if ( !is_u_zero && !is_z_zero )
    dz_u = u_norm * sqrt_pert
    quot = dz_u/z_norm
    @.. z = uprev + quot*z
  elseif !is_u_zero
    dz_u = u_norm * sqrt_pert
    @.. z = uprev + uprev*dz_u
  elseif !is_z_zero
    dz_u = pert
    quot = dz_u/z_norm
    @.. z *= quot
  else
    dz_u = pert
    @.. z = dz_u
  end # endif
  # Start power iteration
  integrator.eigen_est = 0
  for iter in 1:maxiter
    if integrator.alg isa IRKC
      f.f2(fz, z, p, t)
      @.. atmp = fz - cache.du₂
    else
      f(fz, z, p, t)
      integrator.destats.nf += 1
      @.. atmp = fz - fsalfirst
    end
    Δ  = integrator.opts.internalnorm(atmp,t)
    eig_prev = integrator.eigen_est
    integrator.eigen_est = Δ/dz_u * safe
    # Convergence
    if typeof(integrator.alg) <: RKCAlgs # To match the constants given in the paper
      if iter >= 2 && abs(eig_prev - integrator.eigen_est) < max(integrator.eigen_est,1.0/integrator.opts.dtmax)*0.01
        integrator.eigen_est *= 1.2
        # Store the eigenvector
        @.. ccache.zprev = z - uprev
        return true
      end
    else
      if iter >= 2 && abs(eig_prev - integrator.eigen_est) < integrator.eigen_est*0.05
        # Store the eigenvector
        @.. ccache.zprev = z
        return true
      end
    end
    # Next `z`
    if Δ != zero(Δ)
      quot = dz_u/Δ
      @.. z = uprev + quot*atmp
    else
      # An arbitrary change on `z`
      nind = length(uprev)
      if (nind != 1)
        ind = 1 + iter % nind
        z[ind] = uprev[ind] - (z[ind] - uprev[ind])
      else
        z = -z
      end
    end
  end
  return false
end
"""
    choosedeg!(cache) -> nothing

Calculate `ms[mdeg]` (the degree of the Chebyshev polynomial)
and `cache.recind` (the index of recurrence parameters for that
degree), where `recf[recind:(recind+ms[mdeg]-2)]` are the `μ,κ` pairs
for the `mdeg` degree method.
  """
function choosedeg!(cache::T) where T
  isconst = T <: OrdinaryDiffEqConstantCache
  isconst || ( cache = cache.constantcache )
  recind = 0
  @inbounds for i in 1:size(cache.ms,1)
    recind += cache.ms[i]
    if cache.ms[i] > cache.mdeg
      cache.mdeg = i
      cache.recind = recind
      break
    end
  end
  return nothing
end


function choosedeg_SERK!(integrator,cache::T) where T
  isconst = T <: OrdinaryDiffEqConstantCache
  isconst || ( cache = cache.constantcache )
  @unpack ms = cache
  start = 1
  @inbounds for i in 1:size(ms,1)
    if ms[i] < cache.mdeg
      start += ms[i]+1
    else
      cache.start = start
      cache.mdeg = ms[i]
      break
    end
  end
  if integrator.alg isa ESERK5
    if cache.mdeg <= 20
      cache.internal_deg = 2
    elseif cache.mdeg <= 50
      cache.internal_deg = 5
    elseif cache.mdeg <= 100
      cache.internal_deg = 10
    elseif cache.mdeg <= 500
      cache.internal_deg = 50
    elseif cache.mdeg <= 1000
      cache.internal_deg = 100
    elseif cache.mdeg <= 2000
      cache.internal_deg = 200
    end
  end if

  if integrator.alg isa SERK2v2
    cache.internal_deg = cache.deg/10
  end if
  return nothing
end
