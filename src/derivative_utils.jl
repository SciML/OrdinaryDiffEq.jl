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
function calc_J!(integrator, cache::OrdinaryDiffEqConstantCache, is_compos)
  @unpack t,dt,uprev,u,f,p = integrator
  if DiffEqBase.has_jac(f)
    J = f.jac(uprev, p, t)
  else
    error("Jacobian wrapper for constant caches not yet implemented") #TODO
  end
  is_compos && (integrator.eigen_est = opnorm(J, Inf))
  return J
end
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
    WOperator(mass_matrix,gamma,J[;cache=nothing,transform=false])

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
operator (must be a `AbstractDiffEqLinearOperator`).

`WOperator` supports lazy `*` and `mul!` operations, the latter utilizing an
internal cache (can be specified in the constructor; default to regular `Vector`).
It supports all of `AbstractDiffEqLinearOperator`'s interface.
"""
mutable struct WOperator{T,
  MType <: Union{UniformScaling,AbstractMatrix},
  GType <: Real,
  JType <: DiffEqBase.AbstractDiffEqLinearOperator{T},
  CType <: AbstractVector
  } <: DiffEqBase.AbstractDiffEqLinearOperator{T}
  mass_matrix::MType
  gamma::GType
  J::JType
  cache::CType
  transform::Bool
  function WOperator(mass_matrix, gamma, J; cache=nothing, transform=false)
    T = eltype(J)
    # Convert mass_matrix, if needed
    if !isa(mass_matrix, Union{AbstractMatrix,UniformScaling})
      mass_matrix = convert(AbstractMatrix, mass_matrix)
    end
    # Construct the cache, default to regular vector
    if cache == nothing
      cache = Vector{T}(undef, size(J, 1))
    end
    new{T,typeof(mass_matrix),typeof(gamma),typeof(J),typeof(cache)}(mass_matrix,gamma,J,cache,transform)
  end
end
set_gamma!(W::WOperator, gamma) = (W.gamma = gamma; W)
DiffEqBase.update_coefficients!(W::WOperator,u,p,t) = (update_coefficients!(W.J,u,p,t); W)
function Base.convert(::Type{AbstractMatrix}, W::WOperator)
  if W.transform
    W.mass_matrix / W.gamma - convert(AbstractMatrix,W.J)
  else
    W.mass_matrix - W.gamma * convert(AbstractMatrix,W.J)
  end
end
function Base.convert(::Type{Number}, W::WOperator)
  if W.transform
    W.mass_matrix / W.gamma - convert(Number,W.J)
  else
    W.mass_matrix - W.gamma * convert(Number,W.J)
  end
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
    mul!(W.cache, W.J, B)
    Y .-= W.cache
  else
    # Compute mass_matrix * B
    if isa(W.mass_matrix, UniformScaling)
      @. Y = W.mass_matrix.λ * B
    else
      mul!(Y, W.mass_matrix, B)
    end
    # Compute J * B
    mul!(W.cache, W.J, B)
    # Subtract result
    axpy!(-W.gamma, W.cache, Y)
  end
end

"""
    lazy_W(f) -> t/f

Predicate for determining what kind of function supports the use of lazy W.
Also used in cache construction.
"""
lazy_W(f) = false # default
lazy_W(f::ODEFunction) = DiffEqBase.has_jac(f) && isa(f.jac_prototype, DiffEqBase.AbstractDiffEqLinearOperator)
# More to come as support for other *DEFunction is added

function calc_W!(integrator, cache::OrdinaryDiffEqMutableCache, dtgamma, repeat_step, W_transform=false)
  @inbounds begin
    @unpack t,dt,uprev,u,f,p = integrator
    @unpack J,W = cache
    mass_matrix = integrator.f.mass_matrix
    is_compos = typeof(integrator.alg) <: CompositeAlgorithm
    alg = unwrap_alg(integrator, true)

    # calculate W
    new_W = true
    if DiffEqBase.has_invW(f)
      # skip calculation of inv(W) if step is repeated
      !repeat_step && W_transform ? f.invW_t(W, uprev, p, dtgamma, t) :
                                    f.invW(W, uprev, p, dtgamma, t) # W == inverse W
      is_compos && calc_J!(integrator, cache, true)

    else
      # skip calculation of J if step is repeated
      if repeat_step || (alg_can_repeat_jac(alg) &&
                         (!integrator.last_stepfail && cache.newton_iters == 1 &&
                          cache.ηold < alg.new_jac_conv_bound))
        new_jac = false
      else
        new_jac = true
        calc_J!(integrator, cache, is_compos)
      end
      # skip calculation of W if step is repeated
      if !repeat_step && (!alg_can_repeat_jac(alg) ||
                          (integrator.iter < 1 || new_jac ||
                           abs(dt - (t-integrator.tprev)) > 100eps(typeof(integrator.t))))
        if lazy_W(f)
          set_gamma!(W, dtgamma)
          # W.transform = W_transform # necessary?
        else # compute W as a dense matrix
          if W_transform
            for j in 1:length(u), i in 1:length(u)
                W[i,j] = mass_matrix[i,j]/dtgamma - J[i,j]
            end
          else
            for j in 1:length(u), i in 1:length(u)
                W[i,j] = mass_matrix[i,j] - dtgamma*J[i,j]
            end
          end
        end
      else
        new_W = false
      end
    end
    return new_W
  end
end

function calc_W!(integrator, cache::OrdinaryDiffEqConstantCache, dtgamma, repeat_step, W_transform=false)
  @unpack t,uprev,f = integrator
  @unpack uf = cache
  mass_matrix = integrator.sol.prob.mass_matrix
  # calculate W
  uf.t = t
  isarray = typeof(uprev) <: AbstractArray
  iscompo = typeof(integrator.alg) <: CompositeAlgorithm
  if !W_transform
    if lazy_W(f)
      W = WOperator(mass_matrix, dtgamma, deepcopy(f.jac_prototype); transform=false)
    else
      if isarray
        J = ForwardDiff.jacobian(uf,uprev)
      else
        J = ForwardDiff.derivative(uf,uprev)
      end
      W = mass_matrix - dtgamma*J
    end
  else
    if lazy_W(f)
      W = WOperator(mass_matrix, dtgamma, deepcopy(f.jac_prototype); transform=true)
    else
      if isarray
        J = ForwardDiff.jacobian(uf,uprev)
      else
        J = ForwardDiff.derivative(uf,uprev)
      end
      W = mass_matrix*inv(dtgamma) - J
    end
  end
  iscompo && (integrator.eigen_est = isarray ? opnorm(J, Inf) : J)
  W
end

function calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, W_transform)
  calc_tderivative!(integrator, cache, dtd1, repeat_step)
  calc_W!(integrator, cache, dtgamma, repeat_step, W_transform)
  return nothing
end
