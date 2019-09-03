## initialize!

@muladd function initialize!(nlsolver::NLSolver{false}, nlcache::NLNewtonConstantCache, integrator)
  @unpack dt = integrator

  nlcache.invγdt = inv(dt * nlsolver.γ)
  nlcache.tstep = integrator.t + nlsolver.c * dt

  nothing
end

@muladd function initialize!(nlsolver::NLSolver{true}, nlcache::NLNewtonCache, integrator)
  @unpack u,uprev,t,dt,opts = integrator
  @unpack weight = nlsolver.cache

  nlcache.invγdt = inv(dt * nlsolver.γ)
  nlcache.tstep = integrator.t + nlsolver.c * dt 
  calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, u,
                       opts.abstol, opts.reltol, opts.internalnorm, t)

  nothing
end

## compute_step!

"""
    compute_step!(nlsolver::NLSolver, nlcache::Union{NLNewtonConstantCache,NLNewtonCache}, integrator)

Compute next iterate of numerically stable modified Newton iteration
that is specialized for implicit methods.

It solves
```math
z = dt⋅f(tmp + γ⋅z, p, t + c⋅dt)
```
by iterating
```math
(I + (dt⋅γ)J) Δᵏ = dt*f(tmp + γ⋅zᵏ, p, t + c⋅dt) - zᵏ
zᵏ⁺¹ = g(zᵏ) = zᵏ - Δᵏ
```
or, by utilizing a transformation,
```math
W Δᵏ = f(tmp + γ⋅zᵏ, p, t + c⋅dt)/γ - zᵏ/(dt⋅γ)
zᵏ⁺¹ = g(zᵏ) = zᵏ - Δᵏ/(dt⋅γ)
```
where `W = M/(dt⋅γ) - J`, `M` is the mass matrix, `dt` is the step size, `γ` and
`c` are constants, `J` is the Jacobian matrix. This transformation occurs since `c*J` is
O(n^2), while `c*M` is usually much sparser. In the most common case, `M=I`, we
have that `c*M` is O(1) for `I isa UniformScaling`.

# References

M.E.Hoseaa and L.F.Shampine, "Analysis and implementation of TR-BDF2",
Applied Numerical Mathematics, Volume 20, Issues 1–2, February 1996, Pages
21-37.
[doi:10.1016/0168-9274(95)00115-8](https://doi.org/10.1016/0168-9274(95)00115-8).

Ernst Hairer and Gerhard Wanner, "Solving Ordinary Differential
Equations II, Springer Series in Computational Mathematics. ISBN
978-3-642-05221-7. Section IV.8.
[doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7).
"""
@muladd function compute_step!(nlsolver::NLSolver{false}, nlcache::NLNewtonConstantCache, integrator)
  @unpack uprev,t,p,dt,opts = integrator
  @unpack z,tmp,γ,cache = nlsolver
  @unpack tstep,W,invγdt = cache

  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  ustep = @. tmp + γ * z
  if mass_matrix === I
    ztmp = (dt .* f(ustep, p, tstep) .- z) .* invγdt
  else
    ztmp = (dt .* f(ustep, p, tstep) .- mass_matrix * z) .* invγdt
  end
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nf += 1
  end

  dz = _reshape(W \ _vec(ztmp), axes(ztmp))
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nsolve += 1
  end

  atmp = calculate_residuals(dz, uprev, ustep, opts.abstol, opts.reltol, opts.internalnorm, t)
  ndz = opts.internalnorm(atmp, t)

  # compute next iterate  
  nlsolver.ztmp = z .- dz

  ndz
end

@muladd function compute_step!(nlsolver::NLSolver{true}, nlcache::NLNewtonCache, integrator)
  @unpack uprev,t,p,dt,opts = integrator
  @unpack z,dz,tmp,ztmp,k,γ,iter,cache = nlsolver
  @unpack ustep,tstep,atmp,W,new_W,invγdt,linsolve,weight = cache

  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  @.. ustep = tmp + γ * z
  f(k, ustep, p, tstep)
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nf += 1
  end

  if mass_matrix === I
    @.. ztmp = (dt * k - z) * invγdt
  else
    mul!(vec(ztmp), mass_matrix, vec(z))
    @.. ztmp = (dt * k - ztmp) * invγdt
  end

  # update W
  if W isa DiffEqBase.AbstractDiffEqLinearOperator
    update_coefficients!(W, ustep, p, tstep)
  end

  linsolve(vec(dz), W, vec(ztmp), iter == 1 && new_W;
           Pl=DiffEqBase.ScaleVector(weight, true),
           Pr=DiffEqBase.ScaleVector(weight, false), tol=integrator.opts.reltol)
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nsolve += 1
  end

  calculate_residuals!(atmp, dz, uprev, ustep, opts.abstol, opts.reltol, opts.internalnorm, t)
  ndz = opts.internalnorm(atmp, t)

  # compute next iterate
  @.. ztmp = z - dz

  ndz
end