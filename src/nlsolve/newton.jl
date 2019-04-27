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
W Δᵏ = dt⋅f(tmp + γ⋅zᵏ, p, t + c⋅h) - zᵏ
zᵏ⁺¹ = zᵏ + Δᵏ
```

where `W = M - dt⋅γJ`, `M` is the mass matrix, `dt` is the step size, `γ` is a
constant, `J` is the Jacobian matrix.

It returns the tuple `z`, where `z` is the solution.

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
  tstep = t + c*dt
  η = max(nlsolver.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)

  # Newton iteration
  local ndz
  fail_convergence = true
  iter = 0
  while iter < max_iter
    iter += 1
    integrator.destats.nnonliniter += 1

    # evaluate function
    u = @.. tmp + γ * z
    if mass_matrix === I
      ztmp = dt .* f(u, p, tstep) .- z
    else
      ztmp = dt .* f(u, p, tstep) .- mass_matrix * z
    end
    integrator.destats.nf += 1
    if DiffEqBase.has_invW(f)
      dz = _reshape(W * -_vec(ztmp), axes(ztmp)) # Here W is actually invW
    else
      dz = _reshape(W \ _vec(ztmp), axes(ztmp))
    end
    integrator.destats.nsolve += 1

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
    z = z .- dz

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndz < κ && (iter > 1 || iszero(ndz) || !iszero(integrator.success_iter))
      # Newton method converges
      nlsolver.status = η < nlsolver.fast_convergence_cutoff ? FastConvergence : Convergence
      fail_convergence = false
      break
    end
  end

  if fail_convergence
    integrator.destats.nnonlinconvfail += 1
  end
  integrator.force_stepfail = fail_convergence
  nlsolver.ηold = η
  nlsolver.nl_iters = iter
  return z
end

@muladd function nlsolve!(nlsolver::NLSolver, nlcache::NLNewtonCache, integrator)
  @unpack t,dt,uprev,u,p,cache = integrator
  @unpack z,dz,tmp,ztmp,k,κ,c,γ,max_iter = nlsolver
  @unpack W, new_W, W_dt = nlcache
  cache = unwrap_cache(integrator, true)

  # precalculations
  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)
  vecztmp = vec(ztmp); vecz = vec(z); vecdz = vec(dz)
  tstep = t + c*dt
  η = max(nlsolver.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)

  # Newton iteration
  local ndz
  fail_convergence = true
  iter = 0
  while iter < max_iter
    iter += 1
    integrator.destats.nnonliniter += 1

    # evaluate function
    @.. u = tmp + γ*z
    f(k, u, p, tstep)
    integrator.destats.nf += 1
    if mass_matrix === I
      @.. ztmp = dt*k - z
    else
      mul!(vecztmp,mass_matrix,vecz)
      @.. ztmp = dt*k - ztmp
    end
    if DiffEqBase.has_invW(f)
      mul!(vecdz,W,vecztmp) # Here W is actually invW
      @.. vecdz = -vecdz
    else
      cache.linsolve(vecdz,W,vecztmp,iter == 1 && new_W)
    end
    integrator.destats.nsolve += 1

    # compute norm of residuals
    iter > 1 && (ndzprev = ndz)
    #W_dt != dt && (rmul!(dz, 2/(1 + dt / W_dt))) # relaxation
    calculate_residuals!(ztmp, dz, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
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

  if fail_convergence
    integrator.destats.nnonlinconvfail += 1
  end
  integrator.force_stepfail = fail_convergence
  nlsolver.ηold = η
  nlsolver.nl_iters = iter
  return z
end
