"""
    nlsolve!(nlsolver::NLSolver, nlcache::AbstractNLSolverCache, integrator)

Solve
```math
dt⋅f(tmp + γ⋅z, p, t + c⋅dt) = z
```
where `dt` is the step size and `γ` and `c` are constants, and return the solution `z`.
"""
function nlsolve!(nlsolver::NLSolver, nlcache::AbstractNLSolverCache, integrator)
  @unpack max_iter, κ, fast_convergence_cutoff = nlsolver

  initialize!(nlsolver, nlcache, integrator)
  η = initial_η(nlsolver, nlcache, integrator)

  local ndz
  fail_convergence = true
  iter = 0
  while iter < max_iter
    iter += 1
    nlsolver.iter = iter
    if DiffEqBase.has_destats(integrator)
      integrator.destats.nnonliniter += 1
    end

    # compute next step and calculate norm of residuals
    iter > 1 && (ndzprev = ndz)
    ndz = compute_step!(nlsolver, nlcache, integrator)

    # check divergence (not in initial step)
    if iter > 1
      θ = ndz / ndzprev
      ( diverge = θ > 1 ) && ( nlsolver.status = Divergence )
      ( veryslowconvergence = ndz * θ^(max_iter - iter) > κ * (1 - θ) ) && ( nlsolver.status = VerySlowConvergence )
      if diverge || veryslowconvergence
        break
      end
    end

    apply_step!(nlsolver, nlcache, integrator)

    # check for convergence
    iter > 1 && (η = θ / (1 - θ))
    if η * ndz < κ && (iter > 1 || iszero(ndz) || (isnewton(nlsolver) && !iszero(integrator.success_iter)))
      # Newton method converges
      nlsolver.status = η < fast_convergence_cutoff ? FastConvergence : Convergence
      fail_convergence = false
      break
    end

    loopfooter!(nlsolver, nlcache, integrator)
  end

  if fail_convergence && DiffEqBase.has_destats(integrator)
    integrator.destats.nnonlinconvfail += 1
  end
  integrator.force_stepfail = fail_convergence
  nlsolver.ηold = η
  nlsolver.z
end

## default implementations

initialize!(::NLSolver, ::AbstractNLSolverCache, integrator) = nothing

initial_η(nlsolver::NLSolver, ::AbstractNLSolverCache, integrator) = 
  max(nlsolver.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)

function apply_step!(nlsolver::NLSolver{iip}, ::AbstractNLSolverCache, integrator) where iip
  if iip
    @.. nlsolver.z = nlsolver.ztmp
  else
    nlsolver.z = nlsolver.ztmp
  end

  nothing
end

loopfooter!(::NLSolver, ::AbstractNLSolverCache, integrator) = nothing