"""
    nlsolve!(nlsolver::AbstractNLSolver, integrator)

Solve
```math
dt⋅f(tmp + γ⋅z, p, t + c⋅dt) = z
```
where `dt` is the step size and `γ` and `c` are constants, and return the solution `z`.
"""
function nlsolve!(nlsolver::AbstractNLSolver, integrator, cache, γdt, repeat_step)
  @label REDO
  update_W!(integrator, cache, γdt, repeat_step)

  @unpack maxiters, κ, fast_convergence_cutoff = nlsolver

  initialize!(nlsolver, integrator)
  nlsolver.status = Divergence
  θ = nlsolver.ηold

  local ndz
  for iter in 1:maxiters
    nlsolver.iter = iter

    # compute next step and calculate norm of residuals
    iter > 1 && (ndzprev = ndz)
    ndz = compute_step!(nlsolver, integrator)

    # convergence rate
    iter > 1 && (θ = max(0.3 * θ, ndz / ndzprev))

    # check for convergence
    eest = ndz * min(one(θ), θ) / 0.01 # TODO: tune this
    if eest <= one(eest)
      nlsolver.status = Convergence
      break
    end

    # check divergence (not in initial step)
    if iter > 1 && θ > 2
      nlsolver.status = Divergence
      break
    end

    apply_step!(nlsolver, integrator)
  end

  if nlsolver.status == Divergence && !isjacobiancurrent(nlsolver, integrator)
    nlsolver.status = TryAgain
    @goto REDO
  end

  nlsolver.ηold = θ
  postamble!(nlsolver, integrator)
end

## default implementations

initialize!(::AbstractNLSolver, integrator) = nothing

#initial_η(nlsolver::NLSolver, integrator) =
#  max(nlsolver.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)

function apply_step!(nlsolver::NLSolver{algType,iip}, integrator) where {algType,iip}
  if iip
    @.. nlsolver.z = nlsolver.ztmp
  else
    nlsolver.z = nlsolver.ztmp
  end

  nothing
end

function postamble!(nlsolver::NLSolver, integrator)
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nnonliniter += nlsolver.iter

    if nlsolvefail(nlsolver)
      integrator.destats.nnonlinconvfail += 1
    end
  end
  integrator.force_stepfail = nlsolvefail(nlsolver)

  nlsolver.z
end
