## initialize!

@muladd function initialize!(nlsolver::NLSolver,
                             nlcache::Union{NLFunctionalConstantCache,NLFunctionalCache},
                             integrator)
  nlcache.tstep = integrator.t + nlsolver.c * integrator.dt

  nothing
end

@muladd function initialize!(nlsolver::NLSolver,
                             nlcache::Union{NLAndersonConstantCache,NLAndersonCache},
                             integrator)
  nlcache.history = 0
  nlcache.tstep = integrator.t + nlsolver.c * integrator.dt

  nothing
end

## initial_η

function initial_η(nlsolver::NLSolver,
                   nlcache::Union{NLFunctionalCache,NLAndersonCache,
                                  NLFunctionalConstantCache,NLAndersonConstantCache},
                   integrator)
  nlsolver.ηold
end

## compute_step!

function loopfooter!(nlsolver::NLSolver{false}, nlcache::NLAndersonConstantCache, integrator)
  @unpack nl_iters = nlsolver
  @unpack aa_start = nlcache

  # perform Anderson acceleration
  if nl_iters == aa_start
    # update cached values for next step of Anderson acceleration
    nlcache.dzold = nlsolver.dz
    nlcache.z₊old = nlsolver.z
  elseif aa_start < nl_iters < nlsolver.max_iter
    @unpack z,dz = nlsolver
    @unpack Δz₊s,z₊old,dzold,R,Q,γs,history,droptol = nlcache
    # increase size of history
    history += 1

    # remove oldest history if maximum size is exceeded
    max_history = length(Δz₊s)
    if history > max_history
      # circularly shift differences of z₊
      for i in 1:(max_history-1)
        Δz₊s[i] = Δz₊s[i + 1]
      end

      # delete left-most column of QR decomposition
      qrdelete!(Q, R, max_history)

      # update size of history
      history = max_history
    end

    # update history of differences of z₊
    Δz₊s[history] = @.. z - z₊old

    # replace/add difference of residuals as right-most column to QR decomposition
    qradd!(Q, R, _vec(dz .- dzold), history)

    # update cached values
    nlcache.dzold = dz
    nlcache.z₊old = z

    # define current Q and R matrices
    Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

    # check condition (TODO: incremental estimation)
    if droptol !== nothing
      while cond(R) > droptol && history > 1
        qrdelete!(Q, R, history)
        history -= 1
        Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
      end
    end

    # solve least squares problem
    γscur = view(γs, 1:history)
    ldiv!(Rcur, mul!(γscur, Qcur', _vec(dz)))
    if DiffEqBase.has_destats(integrator)
      integrator.destats.nsolve += 1
    end

    # update next iterate
    for i in 1:history
      z = @.. z - γs[i] * Δz₊s[i]
    end
    nlsolver.z = z

    # save updated history
    nlcache.history = history
  end

  nothing
end

function loopfooter!(nlsolver::NLSolver{true}, nlcache::NLAndersonCache, integrator)
  @unpack nl_iters = nlsolver
  @unpack aa_start = nlcache

  # perform Anderson acceleration
  if nl_iters == aa_start
    # update cached values for next step of Anderson acceleration
    @.. nlcache.dzold = nlsolver.dz
    @.. nlcache.z₊old = nlsolver.z
  elseif aa_start < nl_iters < nlsolver.max_iter
    @unpack z,dz = nlsolver
    @unpack z₊old,dzold,Δz₊s,γs,R,Q,history,droptol = nlcache

    # increase size of history
    history += 1

    # remove oldest history if maximum size is exceeded
    max_history = length(Δz₊s)
    if history > max_history
      # circularly shift differences of z₊
      ptr = Δz₊s[1]
      for i in 1:(max_history-1)
        Δz₊s[i] = Δz₊s[i + 1]
      end
      Δz₊s[max_history] = ptr

      # delete left-most column of QR decomposition
      qrdelete!(Q, R, max_history)

      # update size of history
      history = max_history
    end

    # update history of differences of z₊
    @.. Δz₊s[history] = z - z₊old

    # replace/add difference of residuals as right-most column to QR decomposition
    @.. dzold = dz - dzold
    qradd!(Q, R, vec(dzold), history)

    # update cached values
    @.. dzold = dz
    @.. z₊old = z

    # define current Q and R matrices
    Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

    # check condition (TODO: incremental estimation)
    if droptol !== nothing
      while cond(R) > droptol && history > 1
        qrdelete!(Q, R, history)
        history -= 1
        Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
      end
    end

    # solve least squares problem
    γscur = view(γs, 1:history)
    ldiv!(Rcur, mul!(γscur, Qcur', vec(dz)))
    if DiffEqBase.has_destats(integrator)
      integrator.destats.nsolve += 1
    end

    # update next iterate
    for i in 1:history
      @.. z = z - γs[i] * Δz₊s[i]
    end

    # save updated history
    nlcache.history = history
  end

  nothing
end

"""
    compute_step!(nlsolver::NLSolver,
                  nlcache::Union{NLFunctionalCache,NLAndersonCache,
                                 NLFunctionalConstantCache,NLAndersonConstantCache},
                  integrator)

Compute the next step of the fixed-point iteration
```math
g(z) = dt⋅f(tmp + γ⋅z, p, t + c⋅dt),
```
and return the norm of ``g(z) - z``.

# References


Ernst Hairer and Gerhard Wanner, "Solving Ordinary Differential
Equations II, Springer Series in Computational Mathematics. ISBN
978-3-642-05221-7. Section IV.8.
[doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7).
"""
@muladd function compute_step!(nlsolver::NLSolver{false},
                               nlcache::Union{NLAndersonConstantCache,
                                              NLFunctionalConstantCache},
                               integrator)
  @unpack uprev,t,p,dt,opts = integrator
  @unpack z,γ,cache = nlsolver
  @unpack tstep = nlcache

  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  ustep = @.. nlsolver.tmp + γ*z
  if mass_matrix == I
    ztmp = dt .* f(ustep, p, tstep)
    dz = ztmp .- z
  else
    ztmp = _reshape(mass_matrix * _vec(z), axes(z))
    dz = dt .* f(ustep, p, tstep) .- ztmp
    ztmp = z .+ dz
  end
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nf += 1
  end

  # compute norm of residuals
  atmp = calculate_residuals(dz, uprev, ustep, opts.abstol, opts.reltol, opts.internalnorm, t)
  ndz = opts.internalnorm(atmp, t)

  # cache results
  nlsolver.ztmp = ztmp
  nlsolver.dz = dz

  ndz
end

@muladd function compute_step!(nlsolver::NLSolver{true},
                               nlcache::Union{NLFunctionalCache,NLAndersonCache},
                               integrator)
  @unpack uprev,t,p,dt,opts = integrator
  @unpack z,dz,tmp,ztmp,k,γ = nlsolver
  @unpack ustep,tstep,atmp = nlcache

  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  @.. ustep = tmp + γ*z
  f(k, ustep, p, tstep)
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nf += 1
  end
  if mass_matrix == I
    @.. ztmp = dt * k
    @.. dz = ztmp - z
  else
    mul!(vec(ztmp), mass_matrix, vec(z))
    @.. dz = dt * k - ztmp
    @.. ztmp = z + dz
  end

  # compute norm of residuals
  calculate_residuals!(atmp, dz, uprev, ustep, opts.abstol, opts.reltol, opts.internalnorm, t)
  ndz = opts.internalnorm(atmp, t)

  ndz
end