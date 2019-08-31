## initialize!

@muladd function initialize!(nlcache::Union{NLFunctionalConstantCache,NLFunctionalCache},
                             nlsolver::NLSolver{<:NLFunctional}, integrator)
  nlcache.tstep = integrator.t + nlsolver.c * integrator.dt

  nothing
end

@muladd function initialize!(nlcache::Union{NLAndersonConstantCache,NLAndersonCache},
                             nlsolver::NLSolver{<:NLAnderson}, integrator)
  nlcache.history = 0
  nlcache.tstep = integrator.t + nlsolver.c * integrator.dt

  nothing
end

## initial_η

initial_η(nlsolver::NLSolver{<:Union{NLFunctional,NLAnderson}}, integrator) = nlsolver.η

## compute_step!

"""
    compute_step!(nlsolver::NLSolver{<:NLFunctional}, integrator)

Compute the next step of the fixed-point iteration
```math
g(z) = dt⋅f(tmp + γ⋅z, p, t + c⋅dt),
```
and return the norm of ``g(z) - z``.
"""
compute_step!(nlsolver::NLSolver{<:NLFunctional}, integrator) =
  _compute_step_fixedpoint!(nlsolver, integrator)
  
function compute_step!(nlsolver::NLSolver{<:NLAnderson,false}, integrator)
  @unpack iter, alg, cache = nlsolver
  @unpack aa_start = alg

  # perform Anderson acceleration
  previter = iter - 1
  if previter == aa_start
    # update cached values for next step of Anderson acceleration
    cache.gzprev = nlsolver.z
    cache.dzprev = cache.dz
  elseif previter > aa_start
    # actually update the next iterate
    nlsolver.z = anderson(nlsolver.z, cache, integrator)
  end

  # compute step
  _compute_step_fixedpoint!(nlsolver, integrator)
end

function compute_step!(nlsolver::NLSolver{<:NLAnderson,true}, integrator)
  @unpack iter, alg, cache = nlsolver
  @unpack aa_start = alg

  # perform Anderson acceleration
  previter = iter - 1
  if previter == aa_start
    # update cached values for next step of Anderson acceleration
    @.. cache.gzprev = nlsolver.z
    @.. cache.dzprev = cache.dz
  elseif previter > aa_start
    # actually update the next iterate
    anderson!(nlsolver.z, cache, integrator)
  end

  # compute step
  _compute_step_fixedpoint!(nlsolver, integrator)
end

@muladd function _compute_step_fixedpoint!(nlsolver::NLSolver{<:Union{NLFunctional,
                                                                      NLAnderson},false},
                                           integrator)
  @unpack uprev,t,p,dt,opts = integrator
  @unpack z,tmp,γ,cache = nlsolver
  @unpack tstep = cache

  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  ustep = @. tmp + γ * z

  if mass_matrix === I
    gz = dt .* f(ustep, p, tstep)
    dz = gz .- z
  else
    mz = _reshape(mass_matrix * _vec(z), axes(z))
    dz = dt .* f(ustep, p, tstep) .- mz
    gz = z .+ dz
  end
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nf += 1
  end
  
  atmp = calculate_residuals(dz, uprev, ustep, opts.abstol, opts.reltol, opts.internalnorm, t)
  ndz = opts.internalnorm(atmp, t)

  # cache dz if Anderson acceleration is performed
  nlsolver.gz = gz
  if nlsolver.alg isa NLAnderson
    cache.dz = dz
  end

  ndz
end

@muladd function _compute_step_fixedpoint!(nlsolver::NLSolver{<:Union{NLFunctional,
                                                                      NLAnderson},true},
                                           integrator)
  @unpack uprev,t,p,dt,opts = integrator
  @unpack z,gz,tmp,γ,cache = nlsolver
  @unpack ustep,dz,tstep,k,atmp = cache

  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  @.. ustep = tmp + γ * z
  f(k, ustep, p, tstep)
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nf += 1
  end

  if mass_matrix === I
    @.. gz = dt * k
    @.. dz = gz - z
  else
    mul!(vec(dz), mass_matrix, vec(z))
    @.. dz = dt * k - dz
    @.. gz = z + dz
  end

  calculate_residuals!(atmp, dz, uprev, ustep, opts.abstol, opts.reltol, opts.internalnorm, t)
  opts.internalnorm(atmp, t)
end

## resize!

function Base.resize!(nlcache::NLFunctionalCache, i::Int)
  resize!(nlcache.ustep, i)
  resize!(nlcache.dz, i)
  resize!(nlcache.k, i)
  resize!(nlcache.atmp, i)
  nothing
end

function Base.resize!(nlcache::Union{NLAndersonCache,NLAndersonConstantCache}, nlsolver::NLSolver{<:NLAnderson},
  integrator, i::Int)
  resize!(nlcache, nlsolver.alg, i)
end

function Base.resize!(nlcache::NLAndersonCache, nlalg::NLAnderson, i::Int)
  @unpack gzprev,Δgzs = nlcache

  resize!(nlcache.ustep, i)
  resize!(nlcache.dz, i)
  resize!(nlcache.k, i)
  resize!(nlcache.atmp, i)
  resize!(gzprev, i)
  resize!(nlcache.dzprev, i)

  # update history of Anderson cache
  max_history_old = length(Δgzs)
  resize_anderson_history!(nlcache, nlalg, i)

  max_history = length(Δgzs)
  if max_history > max_history_old
    for i in (max_history_old + 1):max_history
      Δgzs[i] = zero(gzprev)
    end
  end

  nothing
end

Base.resize!(nlcache::NLAndersonConstantCache, nlalg::NLAnderson, i::Int) =
  resize_anderson_history!(nlcache, nlalg, i)

function resize_anderson_history!(nlcache, nlalg::NLAnderson, i::Int)
  # determine new maximum history
  max_history_old = length(nlcache.Δgzs)
  max_history = min(nlalg.max_history, nlalg.max_iter, i)

  resize!(nlcache.γs, max_history)
  resize!(nlcache.Δgzs, max_history)

  if max_history != max_history_old
    nlcache.Q = typeof(nlcache.Q)(undef, i, max_history)
    nlcache.R = typeof(nlcache.R)(undef, max_history, max_history)
  end

  nothing
end