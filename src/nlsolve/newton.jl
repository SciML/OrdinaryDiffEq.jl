## initialize!

@muladd function initialize!(nlsolver::NLSolver{<:NLNewton,false}, integrator)
  @unpack dt = integrator
  @unpack cache = nlsolver

  cache.invγdt = inv(dt * nlsolver.γ)
  cache.tstep = integrator.t + nlsolver.c * dt

  nothing
end

@muladd function initialize!(nlsolver::NLSolver{<:NLNewton,true}, integrator)
  @unpack u,uprev,t,dt,opts = integrator
  @unpack cache = nlsolver
  @unpack weight = cache

  cache.invγdt = inv(dt * nlsolver.γ)
  cache.tstep = integrator.t + nlsolver.c * dt
  calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, u,
                       opts.abstol, opts.reltol, opts.internalnorm, t)

  nothing
end

## compute_step!

"""
    compute_step!(nlsolver::NLSolver{<:NLNewton}, integrator)

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
@muladd function compute_step!(nlsolver::NLSolver{<:NLNewton,false}, integrator)
  @unpack uprev,t,p,dt,opts = integrator
  @unpack z,tmp,γ,α,cache = nlsolver
  @unpack tstep,W,invγdt = cache

  f = nlsolve_f(integrator)
  isdae = f isa DAEFunction

  if !isdae
    mass_matrix = integrator.f.mass_matrix
  end

  if isdae
    ustep = @. uprev + z
    dustep = @. (tmp + α * z) * invγdt
    ztmp = f(dustep, ustep, p, t)
  else
    ustep = @. tmp + γ * z
    if mass_matrix === I
      ztmp = (dt .* f(ustep, p, tstep) .- z) .* invγdt
    else
      ztmp = (dt .* f(ustep, p, tstep) .- mass_matrix * z) .* invγdt
    end
  end

  if DiffEqBase.has_destats(integrator)
    integrator.destats.nf += 1
  end

  # update W
  if W isa DiffEqBase.AbstractDiffEqLinearOperator
    W = update_coefficients!(W, ustep, p, tstep)
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

@muladd function compute_step!(nlsolver::NLSolver{<:NLNewton,true}, integrator)
  @unpack uprev,t,p,dt,opts = integrator
  @unpack z,tmp,ztmp,γ,α,iter,cache = nlsolver
  @unpack W_γdt,ustep,tstep,k,atmp,dz,W,new_W,invγdt,linsolve,weight = cache

  f = nlsolve_f(integrator)
  isdae = f isa DAEFunction

  if !isdae
    mass_matrix = integrator.f.mass_matrix
  end

  if isdae
    @.. ztmp = (tmp + α * z) * invγdt
    @.. ustep = uprev + z
    f(k, ztmp, ustep, p, tstep)
  else
    @.. ustep = tmp + γ * z
    f(k, ustep, p, tstep)
  end
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nf += 1
  end

  if isdae
    b = vec(k)
  else
    if mass_matrix === I
      @.. ztmp = (dt * k - z) * invγdt
    else
      mul!(vec(ztmp), mass_matrix, vec(z))
      @.. ztmp = (dt * k - ztmp) * invγdt
    end
    b = vec(ztmp)
  end

  # update W
  if W isa DiffEqBase.AbstractDiffEqLinearOperator
    update_coefficients!(W, ustep, p, tstep)
  end

  linsolve(vec(dz), W, b, iter == 1 && new_W;
           Pl=DiffEqBase.ScaleVector(weight, true),
           Pr=DiffEqBase.ScaleVector(weight, false), tol=integrator.opts.reltol)

  if DiffEqBase.has_destats(integrator)
    integrator.destats.nsolve += 1
  end

  # relaxed Newton
  # Diagonally Implicit Runge-Kutta Methods for Ordinary Differential
  # Equations. A Review, by Christopher A. Kennedy and Mark H. Carpenter
  # page 54.
  if isdae
    γdt = α * invγdt
  else
    γdt = γ * dt
  end

  !(W_γdt ≈ γdt) && (rmul!(dz, 2/(1 + γdt / W_γdt)))

  calculate_residuals!(atmp, dz, uprev, ustep, opts.abstol, opts.reltol, opts.internalnorm, t)
  ndz = opts.internalnorm(atmp, t)

  # compute next iterate
  @.. ztmp = z - dz

  ndz
end

## resize!

function Base.resize!(nlcache::NLNewtonCache, ::AbstractNLSolver, integrator, i::Int)
  resize!(nlcache.ustep, i)
  resize!(nlcache.k, i)
  resize!(nlcache.atmp, i)
  resize!(nlcache.dz, i)
  resize!(nlcache.du1, i)
  if nlcache.jac_config !== nothing
    resize_jac_config!(nlcache.jac_config, i)
  end
  resize!(nlcache.weight, i)

  # resize J and W (or rather create new ones of appropriate size and type)
  resize_J_W!(nlcache, integrator, i)

  nothing
end
