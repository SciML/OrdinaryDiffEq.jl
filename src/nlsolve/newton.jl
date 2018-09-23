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
function (S::NLNewton{false})(integrator)
  nlcache = S.cache
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,tmp,W,κ,tol,c,γ,max_iter,min_iter = nlcache
  mass_matrix = integrator.f.mass_matrix
  #alg = unwrap_alg(integrator, true)
  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
  else
    f = integrator.f
  end
  # precalculations
  κtol = κ*tol

  # initial step of NLNewton iteration
  iter = 1
  tstep = t + c*dt
  u = tmp + γ*z
  b = dt*f(u, p, tstep) - z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz
  # fast pass when `f` is linear
  f isa DiffEqBase.AbstractDiffEqLinearOperator && return (z, nlcache.ηold, iter, false)

  η = max(nlcache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # NLNewton iteration
  fail_convergence = false
  while (do_newton || iter < min_iter) && iter < max_iter
    iter += 1
    u = tmp + γ*z
    b = dt*f(u, p, tstep) - z
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(max_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z = z + dz
  end

  if (iter >= max_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return (z, η, iter, true)
  end

  return (z, η, iter, false)
end

function (S::NLNewton{true})(integrator)
  nlcache = S.cache
  cache = integrator.cache
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,dz,tmp,b,W,κ,tol,k,new_W,c,γ,max_iter,min_iter = nlcache
  mass_matrix = integrator.f.mass_matrix
  #alg = unwrap_alg(integrator, true)
  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
  else
    f = integrator.f
  end
  # precalculations
  κtol = κ*tol
  # initial step of NLNewton iteration
  iter = 1
  tstep = t + c*dt
  @. u = tmp + γ*z
  f(k, u, p, tstep)
  if mass_matrix == I
    @. b = dt*k - z
  else
    mul!(vec(b),mass_matrix,vec(z))
    @. b = dt*k - b
  end
  if DiffEqBase.DiffEqBase.has_invW(f)
    mul!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz

  η = max(nlcache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # NLNewton iteration
  fail_convergence = false
  while (do_newton || iter < min_iter) && iter < max_iter
    iter += 1
    @. u = tmp + γ*z
    f(k, u, p, tstep)
    if mass_matrix == I
      @. b = dt*k - z
    else
      mul!(vec(b),mass_matrix,vec(z))
      @. b = dt*k - b
    end
    if DiffEqBase.DiffEqBase.has_invW(f)
      mul!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(max_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z .+= dz
  end

  if (iter >= max_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return (z, η, iter, true)
  end

  return (z, η, iter, false)
end
