abstract type AbstractNLsolveCache end
#abstract type NLsolveConstantCache <: AbstractNLsolveCache end
#abstract type NLsolveMutableCache <: AbstractNLsolveCache end
struct NLsolveConstantCache{uType,J,uToltype,cType,gType} <: AbstractNLsolveCache
  z::uType
  tmp::uType
  W::J
  κ::uToltype
  tol::uToltype
  γ::gType
  c::cType
end
struct NLsolveMutableCache{rateType,uType,J,uToltype,cType,gType} <: AbstractNLsolveCache
  z::uType
  dz::uType
  tmp::uType
  b::uType
  W::J
  κ::uToltype
  tol::uToltype
  γ::gType
  c::cType
  k::rateType
  new_W::Bool
end

"""
  nlsolve_cache(alg, cache::OrdinaryDiffEqConstantCache, z, tmp, W, γ, c, new_W)
  -> NLsolveConstantCache

Return a wrapper `NLsolveConstantCache` for `cache::OrdinaryDiffEqConstantCache`,
so that `diffeq_nlsolve!` does not need to assume fieldnames of `cache`.
"""
function nlsolve_cache(alg::Union{OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                                  OrdinaryDiffEqNewtonAlgorithm},
                       cache::OrdinaryDiffEqConstantCache, z, tmp, W, γ, c, new_W)
  NLsolveConstantCache(z, tmp, W, cache.κ, cache.tol, γ, c)
end

"""
  nlsolve_cache(alg, cache::OrdinaryDiffEqMutableCache, z, tmp, W, γ, c, new_W)
  -> NLsolveMutableCache

Return a wrapper `NLsolveMutableCache` for `cache::OrdinaryDiffEqMutableCache`,
so that `diffeq_nlsolve!` does not need to assume fieldnames of `cache`.
"""
function nlsolve_cache(alg::Union{OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                                  OrdinaryDiffEqNewtonAlgorithm},
                       cache::OrdinaryDiffEqMutableCache, z, tmp, γ, c, new_W)
  NLsolveMutableCache(z, cache.dz, tmp, cache.b,
                      cache.W, cache.κ, cache.tol, γ, c, cache.k, new_W)
end

"""
  diffeq_nlsolve!(integrator, nlcache, cache, ::Type{Val{:newton}}) -> (z, η,
  iter, fail_convergence)

Perform numerically stable modified Newton iteration that is specialized for
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
function diffeq_nlsolve!(integrator,
                         nlcache::NLsolveConstantCache,
                         cache::OrdinaryDiffEqConstantCache,
                         ::Type{Val{:newton}})
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,tmp,W,κ,tol,c,γ = nlcache
  mass_matrix = integrator.sol.prob.mass_matrix
  alg = unwrap_alg(integrator, true)
  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
  else
    f = integrator.f
  end
  # precalculations
  κtol = κ*tol

  # initial step of Newton iteration
  iter = 1
  tstep = t + c*dt
  u = tmp + γ*z
  b = dt*f(u, p, tstep) - z
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z
    b = dt*f(u, p, tstep) - z
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z = z + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return (z, η, iter, true)
  end

  return (z, η, iter, false)
end

function diffeq_nlsolve!(integrator,
                         nlcache::NLsolveMutableCache,
                         cache::OrdinaryDiffEqMutableCache,
                         ::Type{Val{:newton}})
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,dz,tmp,b,W,κ,tol,k,new_W,c,γ = nlcache
  mass_matrix = integrator.sol.prob.mass_matrix
  alg = unwrap_alg(integrator, true)
  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
  else
    f = integrator.f
  end
  # precalculations
  κtol = κ*tol
  # initial step of Newton iteration
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

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
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
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return (z, η, iter, true)
  end

  return (z, η, iter, false)
end

"""
  diffeq_nlsolve!(integrator, nlcache, cache, ::Type{Val{:functional}}) -> (z,
  η, iter, fail_convergence)

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
function diffeq_nlsolve!(integrator,
                         nlcache::NLsolveConstantCache,
                         cache::OrdinaryDiffEqConstantCache,
                         ::Type{Val{:functional}})
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,tmp,κ,tol,c,γ = nlcache
  mass_matrix = integrator.sol.prob.mass_matrix
  alg = unwrap_alg(integrator, true)
  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
  else
    f = integrator.f
  end
  # precalculations
  tol *= 1//(alg_order(alg)*1000)
  κtol = κ*tol

  # initial step of functional iteration
  iter = 1
  tstep = t + c*dt
  u = tmp + γ*z
  z₊ = dt*f(u, p, tstep)
  dz = z₊ - z
  ndz = integrator.opts.internalnorm(dz)
  z = z₊

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_functional = integrator.success_iter == 0 || η*sqrt(ndz) > κtol

  # functional iteration
  fail_convergence = false
  while (do_functional || iter < alg.min_newton_iter) && iter < alg.max_newton_iter # TODO: rename
    iter += 1
    u = tmp + γ*z
    z₊ = dt*f(u, p, tstep)
    dz = z₊ - z
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_functional = (η*sqrt(ndz) > κtol)
    z = z₊
  end

  if (iter >= alg.max_newton_iter && do_functional) || fail_convergence
    integrator.force_stepfail = true
    return (z, η, iter, true)
  end

  return (z, η, iter, false)
end

function diffeq_nlsolve!(integrator,
                         nlcache::NLsolveMutableCache,
                         cache::OrdinaryDiffEqMutableCache,
                         ::Type{Val{:functional}})
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,dz,tmp,κ,tol,b,k,c,γ = nlcache
  z₊ = b
  mass_matrix = integrator.sol.prob.mass_matrix
  alg = unwrap_alg(integrator, true)
  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
  else
    f = integrator.f
  end
  # precalculations
  tol *= 1//(alg_order(alg)*1000)
  κtol = κ*tol
  # initial step of functional iteration
  iter = 1
  tstep = t + c*dt
  @. u = tmp + γ*z
  f(k, u, p, tstep)
  if mass_matrix == I
    @. z₊ = dt*k
  else
    lmul!(dt, k) # TODO: check unit
    mul!(z₊, mass_matrix, k)
  end
  @. dz = z₊ - z
  @. z = z₊
  ndz = integrator.opts.internalnorm(dz)

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_functional = integrator.success_iter == 0 || η*sqrt(ndz) > κtol

  # Functional iteration
  fail_convergence = false
  while (do_functional || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z
    f(k, u, p, tstep)
    if mass_matrix == I
      @. z₊ = dt*k
    else
      lmul!(dt, k) # TODO: check unit
      mul!(z₊, mass_matrix, k)
    end
    @. dz = z₊ - z
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_functional = (η*sqrt(ndz) > κtol)
    @. z = z₊
  end

  if (iter >= alg.max_newton_iter && do_functional) || fail_convergence
    integrator.force_stepfail = true
    return (z, η, iter, true)
  end

  return (z, η, iter, false)
end
