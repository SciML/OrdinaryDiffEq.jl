"""
    nlsolve!(nlsolver::AbstractNLSolver, integrator)

Solve
```math
dt⋅f(tmp + γ⋅z, p, t + c⋅dt) = z
```
where `dt` is the step size and `γ` and `c` are constants, and return the solution `z`.
"""
function nlsolve!(nlsolver::AbstractNLSolver, integrator)
  @unpack maxiters, κ, fast_convergence_cutoff = nlsolver

  initialize!(nlsolver, integrator)
  η = initial_η(nlsolver, integrator)
  nlsolver.status = SlowConvergence

  local ndz
  for iter in 1:maxiters
    nlsolver.iter = iter

    # compute next step and calculate norm of residuals
    iter > 1 && (ndzprev = ndz)
    ndz = compute_step!(nlsolver, integrator)

    # check divergence (not in initial step)
    if iter > 1
      θ = ndz / ndzprev

      # divergence
      if θ > 1
        nlsolver.status = Divergence
        break
      end

      # very slow convergence
      if ndz * θ^(maxiters - iter) > κ * (1 - θ)
        nlsolver.status = VerySlowConvergence
        break
      end
    end

    apply_step!(nlsolver, integrator)

    # check for convergence
    iter > 1 && (η = θ / (1 - θ))
    if η * ndz < κ && (iter > 1 || iszero(ndz) || (isnewton(nlsolver) && !iszero(integrator.success_iter)))
      nlsolver.status = η < fast_convergence_cutoff ? FastConvergence : Convergence
      break
    end
  end

  nlsolver.ηold = η
  postamble!(nlsolver, integrator)
end

## default implementations

initialize!(::AbstractNLSolver, integrator) = nothing

initial_η(nlsolver::NLSolver, integrator) = 
  max(nlsolver.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)

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