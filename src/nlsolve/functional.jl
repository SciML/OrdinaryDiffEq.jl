"""
  (S::NLFunctional)(integrator) -> (z, η, iter, fail_convergence)

Perform functional iteration that is used by implicit methods, where `z` is the
solution, `η` is used to measure the iteration error (see [^HW96]), `iter` is
the number of iteration, and `fail_convergence` reports whether the algorithm
succeed.  It solves

```math
G(z) = dt⋅f(tmp + γ⋅z, p, t+c⋅h)
z = G(z)
```

by iterating

```math
zᵏ⁺¹ = G(zᵏ).
```

[^HW96]: Ernst Hairer and Gerhard Wanner, "Solving Ordinary Differential
Equations II, Springer Series in Computational Mathematics. ISBN
978-3-642-05221-7. Section IV.8.
[doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7)
"""
function (S::NLFunctional{false})(integrator)
  nlcache = S.cache
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,tmp,κ,tol,c,γ,min_iter,max_iter = nlcache
  mass_matrix = integrator.f.mass_matrix
  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
  else
    f = integrator.f
  end
  # precalculations
  κtol = κ*tol

  # initial step of functional iteration
  iter = 1
  tstep = t + c*dt
  u = tmp + γ*z
  z₊ = dt*f(u, p, tstep)
  dz = z₊ - z
  ndz = integrator.opts.internalnorm(dz)
  z = z₊

  η = nlcache.ηold
  do_functional = true # TODO: this makes `min_iter` ≥ 2

  # functional iteration
  fail_convergence = false
  while (do_functional || iter < min_iter) && iter < max_iter
    iter += 1
    u = tmp + γ*z
    z₊ = dt*f(u, p, tstep)
    dz = z₊ - z
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(max_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_functional = (η*ndz > κtol)
    z = z₊
  end

  if (iter >= max_iter && do_functional) || fail_convergence
    integrator.force_stepfail = true
    return (z, η, iter, true)
  end

  return (z, η, iter, false)
end

function (S::NLFunctional{true})(integrator)
  nlcache = S.cache
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,z₊,b,dz,tmp,κ,tol,k,c,γ,min_iter,max_iter = nlcache
  ztmp = b
  mass_matrix = integrator.f.mass_matrix
  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
  else
    f = integrator.f
  end
  # precalculations
  κtol = κ*tol
  # initial step of functional iteration
  iter = 1
  tstep = t + c*dt
  @. u = tmp + γ*z
  f(k, u, p, tstep)
  if mass_matrix == I
    @. z₊ = dt*k
  else
    @. ztmp = dt*k
    mul!(z₊, mass_matrix, ztmp)
  end
  @. dz = z₊ - z
  @. z = z₊
  ndz = integrator.opts.internalnorm(dz)

  η = nlcache.ηold
  do_functional = true # TODO: this makes `min_iter` ≥ 2

  # NLFunctional iteration
  fail_convergence = false
  while (do_functional || iter < min_iter) && iter < max_iter
    iter += 1
    @. u = tmp + γ*z
    f(k, u, p, tstep)
    if mass_matrix == I
      @. z₊ = dt*k
    else
      @. ztmp = dt*k
      mul!(z₊, mass_matrix, ztmp)
    end
    @. dz = z₊ - z
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(max_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_functional = (η*ndz > κtol)
    @. z = z₊
  end

  if (iter >= max_iter && do_functional) || fail_convergence
    integrator.force_stepfail = true
    return (z, η, iter, true)
  end

  return (z, η, iter, false)
end
