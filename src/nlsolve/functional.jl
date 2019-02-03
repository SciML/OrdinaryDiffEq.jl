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
@muladd function (S::NLFunctional{false})(integrator)
  nlcache = S.cache
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,tmp,κ,tol,c,γ,min_iter,max_iter = nlcache
  mass_matrix = integrator.f.mass_matrix
  f = integrator.f isa SplitFunction ? integrator.f.f1 : integrator.f
  # precalculations
  κtol = κ*tol

  # initial step of functional iteration
  iter = 1
  tstep = t + c*dt
  u = @. tmp + γ * z
  if mass_matrix == I
    z₊ = dt .* f(u, p, tstep)
  else
    mz = mass_matrix * z
    z₊ = dt .* f(u, p, tstep) .- mz .+ z
  end

  # compute initial values for early stopping criterion
  ndz = integrator.opts.internalnorm(z₊ .- z)

  # update solution
  z = z₊

  # check stopping criterion for initial step
  η = nlcache.ηold
  do_functional = true # TODO: this makes `min_iter` ≥ 2

  # functional iteration
  while do_functional && iter < max_iter
    # compute next iterate
    iter += 1
    u = @. tmp + γ * z
    if mass_matrix == I
      z₊ = dt .* f(u, p, tstep)
    else
      mz = mass_matrix * z
      z₊ = dt .* f(u, p, tstep) .- mz .+ z
    end

    # check early stopping criterion
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(z₊ .- z)
    θ = ndz/ndzprev
    if θ ≥ 1 || ndz * θ^(max_iter - iter) > κtol * (1 - θ)
      break
    end

    # update solution
    z = z₊

    # check stopping criterion
    η = θ / (1 - θ) # calculated for possible early stopping
    do_functional = iter < min_iter || η * ndz > κtol
  end

  integrator.force_stepfail = do_functional
  z, η, iter, do_functional
end

@muladd function (S::NLFunctional{true})(integrator)
  nlcache = S.cache
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,z₊,b,dz,tmp,κ,tol,k,c,γ,min_iter,max_iter = nlcache
  mass_matrix = integrator.f.mass_matrix
  f = integrator.f isa SplitFunction ? integrator.f.f1 : integrator.f
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
    @. z₊ = dt*k + z
    mul!(b, mass_matrix, z)
    @. z₊ -= b
  end

  # compute initial values for early stopping criterion
  @. dz = z₊ - z
  ndz = integrator.opts.internalnorm(dz)

  # update solution
  @. z = z₊

  # check stopping criterion for initial step
  η = nlcache.ηold
  do_functional = true # TODO: this makes `min_iter` ≥ 2

  # functional iteration
  while do_functional && iter < max_iter
    # compute next iterate
    iter += 1
    @. u = tmp + γ*z
    f(k, u, p, tstep)
    if mass_matrix == I
      @. z₊ = dt*k
    else
      @. z₊ = dt*k + z
      mul!(b, mass_matrix, z)
      @. z₊ -= b
    end
    @. dz = z₊ - z
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)

    # check early stopping criterion
    θ = ndz/ndzprev
    if θ ≥ 1 || ndz * θ^(max_iter - iter) > κtol * (1 - θ)
      break
    end

    # update solution
    @. z = z₊

    # check stopping criterion
    η = θ / (1 - θ) # calculated for possible early stopping
    do_functional = iter < min_iter || η * ndz > κtol
  end

  integrator.force_stepfail = do_functional
  z, η, iter, do_functional
end
