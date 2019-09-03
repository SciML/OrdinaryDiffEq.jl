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

  local ndz
  fail_convergence = true
  iter = 0
  while iter < maxiters
    iter += 1
    nlsolver.iter = iter
    if DiffEqBase.has_destats(integrator)
      integrator.destats.nnonliniter += 1
    end

    # compute next step and calculate norm of residuals
    iter > 1 && (ndzprev = ndz)
    ndz = compute_step!(nlsolver, integrator)

    # check divergence (not in initial step)
    if iter > 1
      θ = ndz / ndzprev
      ( diverge = θ > 1 ) && ( nlsolver.status = Divergence )
      ( veryslowconvergence = ndz * θ^(maxiters - iter) > κ * (1 - θ) ) && ( nlsolver.status = VerySlowConvergence )
      if diverge || veryslowconvergence
        break
      end
    end

    apply_step!(nlsolver, integrator)

    # check for convergence
    iter > 1 && (η = θ / (1 - θ))
    if η * ndz < κ && (iter > 1 || iszero(ndz) || (isnewton(nlsolver) && !iszero(integrator.success_iter)))
      # Newton method converges
      nlsolver.status = η < fast_convergence_cutoff ? FastConvergence : Convergence
      fail_convergence = false
      break
    end
  end

  if fail_convergence && DiffEqBase.has_destats(integrator)
    integrator.destats.nnonlinconvfail += 1
  end
  integrator.force_stepfail = fail_convergence
  nlsolver.ηold = η
  nlsolver.z
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