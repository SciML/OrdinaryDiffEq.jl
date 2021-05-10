"""
    nlsolve!(nlsolver::AbstractNLSolver, integrator)

Solve
```math
dt⋅f(innertmp + γ⋅z, p, t + c⋅dt) + outertmp = z
```
where `dt` is the step size and `γ` and `c` are constants, and return the solution `z`.
"""
function nlsolve!(nlsolver::AbstractNLSolver, integrator, cache=nothing, repeat_step=false)
  @label REDO
  if isnewton(nlsolver)
    cache === nothing && throw(ArgumentError("cache is not passed to `nlsolve!` when using NLNewton"))
    if nlsolver.method === DIRK
      γW = nlsolver.γ * integrator.dt
    else
      γW = nlsolver.γ * integrator.dt / nlsolver.α
    end
    update_W!(nlsolver, integrator, cache, γW, repeat_step)
  end

  @unpack maxiters, κ, fast_convergence_cutoff = nlsolver

  initialize!(nlsolver, integrator)
  nlsolver.status = Divergence
  η = get_new_W!(nlsolver) ? initial_η(nlsolver, integrator) : nlsolver.ηold

  local ndz
  for iter in 1:maxiters
    nlsolver.iter = iter

    # compute next step and calculate norm of residuals
    iter > 1 && (ndzprev = ndz)
    ndz = compute_step!(nlsolver, integrator)
    if !isfinite(ndz)
      nlsolver.status = Divergence
      break
    end

    # check divergence (not in initial step)
    if iter > 1
      θ = ndz / ndzprev

      # divergence
      if θ > 2
        @warn "newton diverge"
        nlsolver.status = Divergence
        break
      end
    end

    apply_step!(nlsolver, integrator)

    # check for convergence
    iter > 1 && (η = θ / (1 - θ))
    if (iter == 1 && ndz < 1e-5) || (iter > 1 && (η >= zero(η) && η * ndz < κ))
      nlsolver.status = Convergence
      #@info "converge"
      break
    end
  end

  if isnewton(nlsolver) && nlsolver.status == Divergence && !isJcurrent(nlsolver, integrator)
    nlsolver.status = TryAgain
    @goto REDO
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
  setfirststage!(nlsolver, false)
  isnewton(nlsolver) && (nlsolver.cache.firstcall = false)

  nlsolver.z
end
