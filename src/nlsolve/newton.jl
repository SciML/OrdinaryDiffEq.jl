"""
  (S::NLNewton)(integrator) -> (z, η, iter, fail_convergence)

Perform numerically stable modified NLNewton iteration that is specialized for
implicit methods (see [^HS96] and [^HW96]), where `z` is the solution, `η` is
used to measure the iteration error (see [^HW96]), `iter` is the number of
iteration, and `fail_convergence` reports whether the algorithm succeed. It
solves

```math
G(z) = dt⋅f(tmp + γ⋅z, p, t+c⋅h) - z = 0⃗
```

by iterating

```math
W Δᵏ = dt⋅f(tmp + γ⋅zᵏ, p, t+c⋅h) - zᵏ
zᵏ⁺¹ = zᵏ + Δᵏ
```

where `W=M-dt⋅γJ`, `M` is the mass matrix, `dt` is the step size, `γ` is a
constant, and `J` is the Jacobian matrix.

[^HS96]: M.E.Hoseaa and L.F.Shampine, "Analysis and implementation of TR-BDF2",
Applied Numerical Mathematics, Volume 20, Issues 1–2, February 1996, Pages
21-37.
[doi:10.1016/0168-9274(95)00115-8](https://doi.org/10.1016/0168-9274(95)00115-8)

[^HW96]: Ernst Hairer and Gerhard Wanner, "Solving Ordinary Differential
Equations II, Springer Series in Computational Mathematics. ISBN
978-3-642-05221-7. Section IV.8.
[doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7)
"""
@muladd function (S::NLNewton{false})(integrator)
  nlcache = S.cache
  @unpack t,dt,uprev,u,p = integrator
  @unpack z,tmp,W,κ,tol,c,γ,max_iter,min_iter = nlcache
  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)
  # precalculations
  κtol = κ*tol

  # initial step of Newton iteration
  iter = 1
  tstep = t + c*dt
  u = @. tmp + γ * z
  if mass_matrix == I
    ztmp = dt .* f(u, p, tstep) .- z
  else
    ztmp = dt .* f(u, p, tstep) .- mass_matrix * z
  end
  if DiffEqBase.has_invW(f)
    dz = _reshape(W * _vec(ztmp), axes(ztmp)) # Here W is actually invW
  else
    dz = _reshape(W \ _vec(ztmp), axes(ztmp))
  end
  integrator.destats.nsolve += 1

  # compute initial values for early stopping criterion
  ndz = integrator.opts.internalnorm(dz,t)

  # update solution
  z = z .+ dz

  # check stopping criterion for initial step
  η = max(nlcache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = iszero(integrator.success_iter) || iter < min_iter || η * ndz > κtol

  # Newton iteration
  while do_newton && iter < max_iter
    # compute next iterate
    iter += 1
    u = @. tmp + γ * z
    if mass_matrix == I
      ztmp = dt .* f(u, p, tstep) .- z
    else
      ztmp = dt .* f(u, p, tstep) .- mass_matrix * z
    end
    if DiffEqBase.has_invW(f)
      dz = _reshape(W * _vec(ztmp), axes(ztmp)) # Here W is actually invW
    else
      dz = _reshape(W \ _vec(ztmp), axes(ztmp))
    end
    integrator.destats.nsolve += 1

    # check early stopping criterion
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz,t)
    θ = ndz/ndzprev
    if θ ≥ 1 || ndz * θ^(max_iter - iter) > κtol * (1 - θ)
      break
    end

    # update solution
    z = z .+ dz

    # check stopping criterion
    η = θ / (1 - θ) # calculated for possible early stopping
    do_newton = iter < min_iter || η * ndz > κtol
  end

  integrator.force_stepfail = do_newton
  z, η, iter, do_newton
end

@muladd function (S::NLNewton{true})(integrator)
  nlcache = S.cache
  cache = integrator.cache
  @unpack t,dt,uprev,u,p = integrator
  @unpack z,dz,tmp,ztmp,W,κ,tol,k,new_W,c,γ,max_iter,min_iter = nlcache
  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)
  # precalculations
  κtol = κ*tol
  vecztmp = vec(ztmp); vecz = vec(z); vecdz = vec(dz)

  # initial step of Newton iteration
  iter = 1
  tstep = t + c*dt
  @. u = tmp + γ*z
  f(k, u, p, tstep)
  if mass_matrix == I
    @. ztmp = dt*k - z
  else
    mul!(vecztmp,mass_matrix,vecz)
    @. ztmp = dt*k - ztmp
  end
  if DiffEqBase.has_invW(f)
    mul!(vecdz,W,vecztmp) # Here W is actually invW
  else
    cache.linsolve(vecdz,W,vecztmp,new_W)
  end
  integrator.destats.nsolve += 1

  # compute initial values for early stopping criterion
  ndz = integrator.opts.internalnorm(dz,t)

  # update solution
  z .+= dz

  # check stopping criterion for initial step
  η = max(nlcache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = iszero(integrator.success_iter) || iter < min_iter || η * ndz > κtol

  # Newton iteration
  while do_newton && iter < max_iter
    # compute next iterate
    iter += 1
    @. u = tmp + γ*z
    f(k, u, p, tstep)
    if mass_matrix == I
      @. ztmp = dt*k - z
    else
      mul!(vecztmp,mass_matrix,vecz)
      @. ztmp = dt*k - ztmp
    end
    if DiffEqBase.has_invW(f)
      mul!(vecdz,W,vecztmp) # Here W is actually invW
    else
      cache.linsolve(vecdz,W,vecztmp,false)
    end
    integrator.destats.nsolve += 1

    # check early stopping criterion
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz,t)
    θ = ndz/ndzprev
    if θ ≥ 1 || ndz * θ^(max_iter - iter) > κtol * (1 - θ)
      break
    end

    # update solution
    z .+= dz

    # check stopping criterion
    η = θ / (1 - θ) # calculated for possible early stopping
    do_newton = iter < min_iter || η * ndz > κtol
  end

  integrator.force_stepfail = do_newton
  z, η, iter, do_newton
end
