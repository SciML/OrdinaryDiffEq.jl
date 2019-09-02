"""
    nlsolve!(nlsolver::NLSolver, nlcache::Union{NLNewtonCache,NLNewtonConstantCache}, integrator)

Perform numerically stable modified Newton iteration that is specialized for implicit
methods (see [^HS96] and [^HW96]).

It solves

```math
G(z) = dt⋅f(tmp + γ⋅z, p, t + c⋅h) - z = 0
```

by iterating

```math
(I + (dt⋅γ)J) Δᵏ = dt*f(tmp + γ⋅zᵏ, p, t + c⋅h) - zᵏ
zᵏ⁺¹ = zᵏ + Δᵏ
```

or, by utilizing a transformation,

```math
W Δᵏ = f(tmp + γ⋅zᵏ, p, t + c⋅h)/γ - zᵏ/(dt⋅γ)
zᵏ⁺¹ = zᵏ + Δᵏ/(dt⋅γ)
```

where `W = M/(dt⋅γ) - J`, `M` is the mass matrix, `dt` is the step size, `γ` is
a constant, `J` is the Jacobian matrix. This transformation occurs since `c*J` is
O(n^2), while `c*M` is usually much sparser. In the most common case, `M=I`, we
have that `c*M` is O(1) for `I isa UniformScaling`.

This returns `z`, where `z` is the solution.

[^HS96]: M.E.Hoseaa and L.F.Shampine, "Analysis and implementation of TR-BDF2",
Applied Numerical Mathematics, Volume 20, Issues 1–2, February 1996, Pages
21-37.
[doi:10.1016/0168-9274(95)00115-8](https://doi.org/10.1016/0168-9274(95)00115-8)

[^HW96]: Ernst Hairer and Gerhard Wanner, "Solving Ordinary Differential
Equations II, Springer Series in Computational Mathematics. ISBN
978-3-642-05221-7. Section IV.8.
[doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7)
"""
@muladd function nlsolve!(nlsolver::NLSolver, nlcache::NLNewtonConstantCache, integrator)
  @unpack t,dt,uprev,u,p = integrator
  @unpack z,tmp,κ,c,γ,max_iter = nlsolver
  W = nlcache.W

  # precalculations
  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)
  invγdt = inv(dt*γ)
  tstep = t + c*dt
  η = max(nlsolver.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)

  # Newton iteration
  local ndz
  fail_convergence = true
  iter = 0
  while iter < max_iter
    iter += 1
    if DiffEqBase.has_destats(integrator)
      integrator.destats.nnonliniter += 1
    end

    # evaluate function
    u = @.. tmp + γ * z
    if mass_matrix === I
      ztmp = (dt .* f(u, p, tstep) .- z) .* invγdt
    else
      ztmp = (dt .* f(u, p, tstep) .- mass_matrix * z) .* invγdt
    end
    if DiffEqBase.has_destats(integrator)
      integrator.destats.nf += 1
    end

    dz = _reshape(W \ _vec(ztmp), axes(ztmp))

    if DiffEqBase.has_destats(integrator)
      integrator.destats.nsolve += 1
    end

    # compute norm of residuals
    iter > 1 && (ndzprev = ndz)
    atmp = calculate_residuals(dz, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    ndz = integrator.opts.internalnorm(atmp, t)

    # check divergence (not in initial step)
    if iter > 1
      θ = ndz / ndzprev
      ( diverge = θ > 1 ) && ( nlsolver.status = Divergence )
      ( veryslowconvergence = ndz * θ^(max_iter - iter) > κ * (1 - θ) ) && ( nlsolver.status = VerySlowConvergence )
      if diverge || veryslowconvergence
        # Newton method diverges
        break
      end
    end

    # update solution
    z = @.. z - dz

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndz < κ && (iter > 1 || iszero(ndz) || !iszero(integrator.success_iter))
      # Newton method converges
      nlsolver.status = η < nlsolver.fast_convergence_cutoff ? FastConvergence : Convergence
      fail_convergence = false
      break
    end
  end

  if fail_convergence && DiffEqBase.has_destats(integrator)
    integrator.destats.nnonlinconvfail += 1
  end
  integrator.force_stepfail = fail_convergence
  nlsolver.ηold = η
  nlsolver.nl_iters = iter
  return z
end

@muladd function nlsolve!(nlsolver::NLSolver, nlcache::NLNewtonCache, integrator)
  @unpack t,dt,uprev,u,p,cache = integrator
  @unpack z,dz,tmp,ztmp,k,κ,c,γ,max_iter,weight = nlsolver
  @unpack W, new_W, W_dt = nlcache
  cache = unwrap_cache(integrator, true)
  calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
  lintol = integrator.opts.reltol

  # precalculations
  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)
  invγdt = inv(dt*γ)
  vecztmp = vec(ztmp); vecz = vec(z); vecdz = vec(dz)
  tstep = t + c*dt
  η = max(nlsolver.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)

  # Newton iteration
  local ndz
  fail_convergence = true
  iter = 0
  while iter < max_iter
    iter += 1
    if DiffEqBase.has_destats(integrator)
      integrator.destats.nnonliniter += 1
    end

    # evaluate function
    @.. dz = tmp + γ*z
    f(k, dz, p, tstep)
    if DiffEqBase.has_destats(integrator)
      integrator.destats.nf += 1
    end
    if mass_matrix === I
      @.. ztmp = (dt*k - z) * invγdt
    else
      mul!(vecztmp,mass_matrix,vecz)
      @.. ztmp = (dt*k - ztmp) * invγdt
    end

    if W isa DiffEqBase.AbstractDiffEqLinearOperator
      update_coefficients!(W,dz,p,tstep)
    end
    nlsolver.linsolve(vecdz,W,vecztmp,iter == 1 && new_W; Pl=DiffEqBase.ScaleVector(weight, true), Pr=DiffEqBase.ScaleVector(weight, false), tol=lintol)

    if DiffEqBase.has_destats(integrator)
      integrator.destats.nsolve += 1
    end

    # compute norm of residuals
    iter > 1 && (ndzprev = ndz)
    #W_dt != dt && (rmul!(dz, 2/(1 + dt / W_dt))) # relaxation
    @.. ztmp = tmp + γ*z
    calculate_residuals!(ztmp, dz, uprev, ztmp, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    ndz = integrator.opts.internalnorm(ztmp, t)

    # check divergence (not in initial step)
    if iter > 1
      θ = ndz / ndzprev
      ( diverge = θ > 1 ) && ( nlsolver.status = Divergence )
      ( veryslowconvergence = ndz * θ^(max_iter - iter) > κ * (1 - θ) ) && ( nlsolver.status = VerySlowConvergence )
      if diverge || veryslowconvergence
        break
      end
    end

    # update solution
    @.. z = z - dz

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndz < κ && (iter > 1 || iszero(ndz) || !iszero(integrator.success_iter))
      # Newton method converges
      nlsolver.status = η < nlsolver.fast_convergence_cutoff ? FastConvergence : Convergence
      fail_convergence = false
      break
    end
  end

  if fail_convergence && DiffEqBase.has_destats(integrator)
    integrator.destats.nnonlinconvfail += 1
  end
  integrator.force_stepfail = fail_convergence
  nlsolver.ηold = η
  nlsolver.nl_iters = iter
  return z
end
