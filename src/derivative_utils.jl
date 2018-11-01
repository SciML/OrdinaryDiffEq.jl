function calc_tderivative!(integrator, cache, dtd1, repeat_step)
  @inbounds begin
    @unpack t,dt,uprev,u,f,p = integrator
    @unpack du2,fsalfirst,dT,tf,linsolve_tmp = cache

    # Time derivative
    if !repeat_step # skip calculation if step is repeated
      if DiffEqBase.has_tgrad(f)
        f.tgrad(dT, uprev, p, t)
      else
        tf.uprev = uprev
        tf.p = p
        derivative!(dT, tf, t, du2, integrator, cache.grad_config)
      end
    end

    f(fsalfirst, uprev, p, t)
    @. linsolve_tmp = fsalfirst + dtd1*dT
  end
  nothing
end

function calc_tderivative(integrator, cache)
  @unpack t,dt,uprev,u,f,p = integrator

  # Time derivative
  if DiffEqBase.has_tgrad(f)
    dT = f.tgrad(uprev, p, t)
  else
    tf = cache.tf
    tf.u = uprev
    tf.p = p
    dT = derivative(tf, t, integrator)
  end
  dT
end

"""
    calc_J(integrator,cache,is_compos)

Interface for calculating the jacobian.

For constant caches, a new jacobian object is returned whereas for mutable
caches `cache.J` is updated. In both cases, if `integrator.f` has a custom
jacobian update function, then it will be called for the update. Otherwise,
either ForwardDiff or finite difference will be used depending on the
`jac_config` of the cache.
"""
function calc_J(integrator, cache::OrdinaryDiffEqConstantCache, is_compos)
  @unpack t,dt,uprev,u,f,p = integrator
  if DiffEqBase.has_jac(f)
    J = f.jac(uprev, p, t)
  else
    J = jacobian(cache.uf,uprev,integrator)
  end
  is_compos && (integrator.eigen_est = opnorm(J, Inf))
  return J
end

"""
    calc_J!(integrator,cache,is_compos)

Interface for calculating the jacobian.

For constant caches, a new jacobian object is returned whereas for mutable
caches `cache.J` is updated. In both cases, if `integrator.f` has a custom
jacobian update function, then it will be called for the update. Otherwise,
either ForwardDiff or finite difference will be used depending on the
`jac_config` of the cache.
"""
function calc_J!(integrator, cache::OrdinaryDiffEqMutableCache, is_compos)
  @unpack t,dt,uprev,u,f,p = integrator
  J = cache.J
  if DiffEqBase.has_jac(f)
    f.jac(J, uprev, p, t)
  else
    @unpack du1,uf,jac_config = cache
    uf.t = t
    uf.p = p
    jacobian!(J, uf, uprev, du1, integrator, jac_config)
  end
  is_compos && (integrator.eigen_est = opnorm(J, Inf))
end

"""
    WOperator(mass_matrix,gamma,J[;transform=false])

A linear operator that represents the W matrix of an ODEProblem, defined as

```math
W = MM - \\gamma J
```

or, if `transform=true`:

```math
W = \\frac{1}{\\gamma}MM - J
```

where `MM` is the mass matrix (a regular `AbstractMatrix` or a `UniformScaling`),
`γ` is a real number proportional to the time step, and `J` is the Jacobian
operator (must be a `AbstractDiffEqLinearOperator`). A `WOperator` can also be
constructed using a `*DEFunction` directly as

    WOperator(f,gamma[;transform=false])

`f` needs to have a jacobian and `jac_prototype`, but the prototype does not need
to be a diffeq operator --- it will automatically be converted to one.

`WOperator` supports lazy `*` and `mul!` operations, the latter utilizing an
internal cache (can be specified in the constructor; default to regular `Vector`).
It supports all of `AbstractDiffEqLinearOperator`'s interface.
"""
mutable struct WOperator{T,
  MType <: Union{UniformScaling,AbstractMatrix},
  GType <: Real,
  JType <: DiffEqBase.AbstractDiffEqLinearOperator{T}
  } <: DiffEqBase.AbstractDiffEqLinearOperator{T}
  mass_matrix::MType
  gamma::GType
  J::JType
  transform::Bool       # true => W = mm/gamma - J; false => W = mm - gamma*J
  _func_cache           # cache used in `mul!`
  _concrete_form         # non-lazy form (matrix/number) of the operator
  WOperator(mass_matrix, gamma, J; transform=false) = new{eltype(J),typeof(mass_matrix),
    typeof(gamma),typeof(J)}(mass_matrix,gamma,J,transform,nothing,nothing)
end
function WOperator(f::DiffEqBase.AbstractODEFunction, gamma; transform=false)
  @assert DiffEqBase.has_jac(f) "f needs to have an associated jacobian"
  if isa(f, Union{SplitFunction, DynamicalODEFunction})
    error("WOperator does not support $(typeof(f)) yet")
  end
  # Convert mass matrix, if needed
  mass_matrix = f.mass_matrix
  if !isa(mass_matrix, Union{AbstractMatrix,UniformScaling})
    mass_matrix = convert(AbstractMatrix, mass_matrix)
  end
  # Convert jacobian, if needed
  J = deepcopy(f.jac_prototype)
  if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
    J = DiffEqArrayOperator(J; update_func=f.jac)
  end
  return WOperator(mass_matrix, gamma, J; transform=transform)
end

set_gamma!(W::WOperator, gamma) = (W.gamma = gamma; W)
DiffEqBase.update_coefficients!(W::WOperator,u,p,t) = (update_coefficients!(W.J,u,p,t); W)
function Base.convert(::Type{AbstractMatrix}, W::WOperator)
  if W._concrete_form == nothing
    # Allocating
    if W.transform
      W._concrete_form = W.mass_matrix / W.gamma - convert(AbstractMatrix,W.J)
    else
      W._concrete_form = W.mass_matrix - W.gamma * convert(AbstractMatrix,W.J)
    end
  else
    # Non-allocating
    if W.transform
      rmul!(copyto!(W._concrete_form, W.mass_matrix), 1/W.gamma)
      axpy!(-1, convert(AbstractMatrix,W.J), W._concrete_form)
    else
      copyto!(W._concrete_form, W.mass_matrix)
      axpy!(-W.gamma, convert(AbstractMatrix,W.J), W._concrete_form)
    end
  end
  W._concrete_form
end
function Base.convert(::Type{Number}, W::WOperator)
  if W.transform
    W._concrete_form = W.mass_matrix / W.gamma - convert(Number,W.J)
  else
    W._concrete_form = W.mass_matrix - W.gamma * convert(Number,W.J)
  end
  W._concrete_form
end
Base.size(W::WOperator, args...) = size(W.J, args...)
function Base.getindex(W::WOperator, i::Int)
  if W.transform
    W.mass_matrix[i] / W.gamma - W.J[i]
  else
    W.mass_matrix[i] - W.gamma * W.J[i]
  end
end
function Base.getindex(W::WOperator, I::Vararg{Int,N}) where {N}
  if W.transform
    W.mass_matrix[I...] / W.gamma - W.J[I...]
  else
    W.mass_matrix[I...] - W.gamma * W.J[I...]
  end
end
function Base.:*(W::WOperator, x::Union{AbstractVecOrMat,Number})
  if W.transform
    (W.mass_matrix*x) / W.gamma - W.J*x
  else
    W.mass_matrix*x - W.gamma * (W.J*x)
  end
end
function Base.:\(W::WOperator, x::Union{AbstractVecOrMat,Number})
  if size(W) == () # scalar operator
    convert(Number,W) \ x
  else
    convert(AbstractMatrix,W) \ x
  end
end
function LinearAlgebra.mul!(Y::AbstractVecOrMat, W::WOperator, B::AbstractVecOrMat)
  if W._func_cache == nothing
    # Allocate cache only if needed
    W._func_cache = Vector{eltype(W)}(undef, size(Y, 1))
  end
  if W.transform
    # Compute mass_matrix * B
    if isa(W.mass_matrix, UniformScaling)
      a = W.mass_matrix.λ / W.gamma
      @. Y = a * B
    else
      mul!(Y, W.mass_matrix, B)
      lmul!(1/W.gamma, Y)
    end
    # Compute J * B and subtract
    mul!(W._func_cache, W.J, B)
    Y .-= W._func_cache
  else
    # Compute mass_matrix * B
    if isa(W.mass_matrix, UniformScaling)
      @. Y = W.mass_matrix.λ * B
    else
      mul!(Y, W.mass_matrix, B)
    end
    # Compute J * B
    mul!(W._func_cache, W.J, B)
    # Subtract result
    axpy!(-W.gamma, W._func_cache, Y)
  end
end

function calc_W!(integrator, cache::OrdinaryDiffEqMutableCache, dtgamma, repeat_step, W_transform=false)
  @inbounds begin
    @unpack t,dt,uprev,u,f,p = integrator
    @unpack J,W = cache
    mass_matrix = integrator.f.mass_matrix
    is_compos = typeof(integrator.alg) <: CompositeAlgorithm
    alg = unwrap_alg(integrator, true)
    isnewton = !(typeof(alg) <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm ||
                 typeof(alg) <: OrdinaryDiffEqRosenbrockAlgorithm)
    isnewton && ( nlcache = cache.nlsolve.cache; @unpack ηold,nl_iters = cache.nlsolve.cache)

    # fast pass
    # we only want to factorize the linear operator once
    new_jac = true
    new_W = true
    if (f isa ODEFunction && islinear(f.f)) || (f isa SplitFunction && islinear(f.f1.f))
      new_jac = false
      @goto J2W # Jump to W calculation directly, because we already have J
    end

    # calculate W
    if DiffEqBase.has_invW(f)
      # skip calculation of inv(W) if step is repeated
      !repeat_step && W_transform ? f.invW_t(W, uprev, p, dtgamma, t) :
                                    f.invW(W, uprev, p, dtgamma, t) # W == inverse W
      is_compos && calc_J!(integrator, cache, true)

    elseif DiffEqBase.has_jac(f) && f.jac_prototype != nothing
      # skip calculation of J if step is repeated
      if repeat_step || (alg_can_repeat_jac(alg) &&
                         (!integrator.last_stepfail && nl_iters == 1 &&
                          ηold < alg.new_jac_conv_bound))
        new_jac = false
      else
        new_jac = true
        DiffEqBase.update_coefficients!(W,uprev,p,t)
      end
      @label J2W
      # skip calculation of W if step is repeated
      if !repeat_step && (!alg_can_repeat_jac(alg) ||
                          (integrator.iter < 1 || new_jac ||
                           abs(dt - (t-integrator.tprev)) > 100eps(typeof(integrator.t))))
        W.transform = W_transform
        set_gamma!(W, dtgamma)
      else
        new_W = false
      end
    else # concrete W using jacobian from `calc_J!`
      # skip calculation of J if step is repeated
      if repeat_step || (alg_can_repeat_jac(alg) &&
                         (!integrator.last_stepfail && nl_iters == 1 &&
                          ηold < alg.new_jac_conv_bound))
        new_jac = false
      else
        new_jac = true
        calc_J!(integrator, cache, is_compos)
      end
      # skip calculation of W if step is repeated
      if !repeat_step && (!alg_can_repeat_jac(alg) ||
                          (integrator.iter < 1 || new_jac ||
                           abs(dt - (t-integrator.tprev)) > 100eps(typeof(integrator.t))))
        if W_transform
          for j in 1:length(u), i in 1:length(u)
              W[i,j] = mass_matrix[i,j]/dtgamma - J[i,j]
          end
        else
          for j in 1:length(u), i in 1:length(u)
              W[i,j] = mass_matrix[i,j] - dtgamma*J[i,j]
          end
        end
      else
        new_W = false
      end
    end
    isnewton && ( nlcache.new_W = new_W )
  end
  return nothing
end

function calc_W!(integrator, cache::OrdinaryDiffEqConstantCache, dtgamma, repeat_step, W_transform=false)
  @unpack t,uprev,p,f = integrator
  @unpack uf = cache
  mass_matrix = integrator.f.mass_matrix
  isarray = typeof(uprev) <: AbstractArray
  # calculate W
  uf.t = t
  is_compos = typeof(integrator.alg) <: CompositeAlgorithm
  if (f isa ODEFunction && islinear(f.f)) || (f isa SplitFunction && islinear(f.f1.f))
    J = f.f1.f
    W = WOperator(mass_matrix, dtgamma, J; transform=W_transform)
  elseif DiffEqBase.has_jac(f)
    J = f.jac(uprev, p, t)
    if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
      J = DiffEqArrayOperator(J)
    end
    W = WOperator(mass_matrix, dtgamma, J; transform=W_transform)
  else
    J = calc_J(integrator, cache, is_compos)
    W = W_transform ? mass_matrix*inv(dtgamma) - J :
                      mass_matrix - dtgamma*J
  end
  is_compos && (integrator.eigen_est = isarray ? opnorm(J, Inf) : J)
  W
end

function calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, W_transform)
  calc_tderivative!(integrator, cache, dtd1, repeat_step)
  calc_W!(integrator, cache, dtgamma, repeat_step, W_transform)
  return nothing
end
