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
    integrator.destats.nf += 1
    @.. linsolve_tmp = fsalfirst + dtd1*dT
  end
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
    calc_J(integrator, cache)

Return a new Jacobian object.

If `integrator.f` has a custom Jacobian update function, then it will be called. Otherwise,
either automatic or finite differencing will be used depending on the `uf` object of the
cache.
"""
function calc_J(integrator, cache)
  @unpack t,uprev,f,p,alg = integrator

  if DiffEqBase.has_jac(f)
    J = f.jac(uprev, p, t)
  else
    @unpack uf = cache

    uf.f = nlsolve_f(f, alg)
    uf.p = p
    uf.t = t

    J = jacobian(uf, uprev, integrator)
  end

  integrator.destats.njacs += 1

  if alg isa CompositeAlgorithm
    integrator.eigen_est = opnorm(J, Inf)
  end

  J
end

"""
    calc_J!(J, integrator, cache) -> J

Update the Jacobian object `J`.

If `integrator.f` has a custom Jacobian update function, then it will be called. Otherwise,
either automatic or finite differencing will be used depending on the `cache`.
"""
function calc_J!(J, integrator, cache)
  @unpack t,uprev,f,p,alg = integrator

  if DiffEqBase.has_jac(f)
    f.jac(J, uprev, p, t)
  else
    @unpack du1, uf, jac_config = cache

    uf.f = nlsolve_f(f, alg)
    uf.t = t
    uf.p = p

    jacobian!(J, uf, uprev, du1, integrator, jac_config)
  end

  integrator.destats.njacs += 1

  if alg isa CompositeAlgorithm
    integrator.eigen_est = opnorm(J, Inf)
  end

  return nothing
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
  inplace::Bool
  _func_cache           # cache used in `mul!`
  _concrete_form         # non-lazy form (matrix/number) of the operator
  WOperator(mass_matrix, gamma, J, inplace; transform=false) = new{eltype(J),typeof(mass_matrix),
    typeof(gamma),typeof(J)}(mass_matrix,gamma,J,transform,inplace,nothing,nothing)
end
function WOperator(f::DiffEqBase.AbstractODEFunction, gamma, inplace; transform=false)
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
  return WOperator(mass_matrix, gamma, J, inplace; transform=transform)
end

set_gamma!(W::WOperator, gamma) = (W.gamma = gamma; W)
DiffEqBase.update_coefficients!(W::WOperator,u,p,t) = (update_coefficients!(W.J,u,p,t); W)
function Base.convert(::Type{AbstractMatrix}, W::WOperator)
  if W._concrete_form === nothing || !W.inplace
    # Allocating
    if W.transform
      W._concrete_form = -W.mass_matrix / W.gamma + convert(AbstractMatrix,W.J)
    else
      W._concrete_form = -W.mass_matrix + W.gamma * convert(AbstractMatrix,W.J)
    end
  else
    # Non-allocating
    _W = W._concrete_form
    J = convert(AbstractMatrix,W.J)
    if W.transform
      if _W isa Diagonal # axpby doesn't specialize on Diagonal matrix
        @inbounds for i in axes(W._concrete_form, 1)
          _W[i, i] = J[i, i] - inv(W.gamma) * W.mass_matrix[i, i]
        end
      else
        copyto!(_W, W.mass_matrix)
        axpby!(one(W.gamma), J, -inv(W.gamma), _W)
      end
    else
      if _W isa Diagonal # axpby doesn't specialize on Diagonal matrix
        @inbounds for i in axes(W._concrete_form, 1)
          _W[i, i] = W.gamma*J[i, i] - W.mass_matrix[i, i]
        end
      else
        copyto!(_W, W.mass_matrix)
        axpby!(W.gamma, J, -one(W.gamma), W._concrete_form)
      end
    end
  end
  W._concrete_form
end
function Base.convert(::Type{Number}, W::WOperator)
  if W.transform
    W._concrete_form = -W.mass_matrix / W.gamma + convert(Number,W.J)
  else
    W._concrete_form = -W.mass_matrix + W.gamma * convert(Number,W.J)
  end
  W._concrete_form
end
Base.size(W::WOperator, args...) = size(W.J, args...)
function Base.getindex(W::WOperator, i::Int)
  if W.transform
    -W.mass_matrix[i] / W.gamma + W.J[i]
  else
    -W.mass_matrix[i] + W.gamma * W.J[i]
  end
end
function Base.getindex(W::WOperator, I::Vararg{Int,N}) where {N}
  if W.transform
    -W.mass_matrix[I...] / W.gamma + W.J[I...]
  else
    -W.mass_matrix[I...] + W.gamma * W.J[I...]
  end
end
function Base.:*(W::WOperator, x::Union{AbstractVecOrMat,Number})
  if W.transform
    (W.mass_matrix*x) / -W.gamma + W.J*x
  else
    -W.mass_matrix*x + W.gamma * (W.J*x)
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
  if W._func_cache === nothing
    # Allocate cache only if needed
    W._func_cache = Vector{eltype(B)}(undef, size(Y, 1))
  end
  if W.transform
    # Compute mass_matrix * B
    if isa(W.mass_matrix, UniformScaling)
      a = -W.mass_matrix.λ / W.gamma
      @.. Y = a * B
    else
      mul!(Y, W.mass_matrix, B)
      lmul!(-1/W.gamma, Y)
    end
    # Compute J * B and add
    mul!(W._func_cache, W.J, B)
    Y .+= W._func_cache
  else
    # Compute mass_matrix * B
    if isa(W.mass_matrix, UniformScaling)
      @.. Y = W.mass_matrix.λ * B
    else
      mul!(Y, W.mass_matrix, B)
    end
    # Compute J * B
    mul!(W._func_cache, W.J, B)
    # Add result
    axpby!(W.gamma, W._func_cache, -one(W.gamma), Y)
  end
end

"""
    islinearfunction(integrator) -> Tuple{Bool,Bool}

return the tuple `(is_linear_wrt_odealg, islinearodefunction)`.
"""
islinearfunction(integrator) = islinearfunction(integrator.f, integrator.alg)

"""
    islinearfunction(f, alg) -> Tuple{Bool,Bool}

return the tuple `(is_linear_wrt_odealg, islinearodefunction)`.
"""
function islinearfunction(f, alg)::Tuple{Bool,Bool}
  isode = f isa ODEFunction && islinear(f.f)
  islin = isode || (alg isa SplitAlgorithms && f isa SplitFunction && islinear(f.f1.f))
  return islin, isode
end

#function do_newJ(integrator, alg, cache, repeat_step)::Bool # any changes here need to be reflected in FIRK
#  integrator.iter <= 1 && return true
#  repeat_step && return false
#  first(islinearfunction(integrator)) && return false
#  integrator.opts.adaptive || return true
#  alg_can_repeat_jac(alg) || return true
#  integrator.u_modified && return true
#  # below is Newton specific logic, so we return non-Newton algs here
#  alg isa NewtonAlgorithm || return true
#  isfirk = alg isa RadauIIA5
#  nlstatus = isfirk ? cache.status : get_status(cache.nlsolver)
#  #@show isjacobiancurrent(integrator, cache.nlsolver)
#  nlsolvefail(nlstatus) && return true
#  # no reuse when the cutoff is 0
#  fast_convergence_cutoff = isfirk ? alg.fast_convergence_cutoff : cache.nlsolver.fast_convergence_cutoff
#  iszero(fast_convergence_cutoff) && return true
#  # reuse J when there is fast convergence
#  fastconvergence = nlstatus === FastConvergence
#  return !fastconvergence
#end
#
#function do_newW(integrator, nlsolver, new_jac, W_dt)::Bool # any changes here need to be reflected in FIRK
#  nlsolver === nothing && return true
#  new_jac && return true
#  # reuse W when the change in stepsize is small enough
#  dt = integrator.dt
#  smallstepchange = abs((dt-W_dt)/W_dt) <= get_new_W_dt_cutoff(nlsolver)
#  return !smallstepchange
#end

function do_newJW(integrator, alg, nlsolver, repeat_step)::NTuple{2,Bool}
  repeat_step && return false, false
  (integrator.iter <= 1 && (isdefined(nlsolver, :iter) && nlsolver.iter <= 1)) && true, true
  W_dt = nlsolver.cache.W_dt
  smallstepchange = abs((integrator.dt-W_dt)/W_dt) <= get_new_W_dt_cutoff(nlsolver)
  jbad = nlsolver.status === TryAgain && smallstepchange
  return jbad, !smallstepchange
end

@noinline _throwWJerror(W, J) = throw(DimensionMismatch("W: $(axes(W)), J: $(axes(J))"))
@noinline _throwWMerror(W, mass_matrix) = throw(DimensionMismatch("W: $(axes(W)), mass matrix: $(axes(mass_matrix))"))
@noinline _throwJMerror(J, mass_matrix) = throw(DimensionMismatch("J: $(axes(J)), mass matrix: $(axes(mass_matrix))"))

@inline function jacobian2W!(W::AbstractMatrix, mass_matrix::MT, dtgamma::Number, J::AbstractMatrix, W_transform::Bool)::Nothing where MT
  # check size and dimension
  iijj = axes(W)
  @boundscheck (iijj === axes(J) && length(iijj) === 2) || _throwWJerror(W, J)
  mass_matrix isa UniformScaling || @boundscheck axes(mass_matrix) === axes(W) || _throwWMerror(W, mass_matrix)
  @inbounds if W_transform
    invdtgamma = inv(dtgamma)
    if MT <: UniformScaling
      copyto!(W, J)
      idxs = diagind(W)
      λ = -mass_matrix.λ
      @.. @view(W[idxs]) = muladd(λ, invdtgamma, @view(J[idxs]))
    else
      @.. W = muladd(-mass_matrix, invdtgamma, J)
    end
  else
    if MT <: UniformScaling
      idxs = diagind(W)
      @.. W = dtgamma*J
      λ = -mass_matrix.λ
      @.. @view(W[idxs]) = @view(W[idxs]) + λ
    else
      @.. W = muladd(dtgamma, J, -mass_matrix)
    end
  end
  return nothing
end

@inline function jacobian2W(mass_matrix::MT, dtgamma::Number, J::AbstractMatrix, W_transform::Bool)::Nothing where MT
  # check size and dimension
  mass_matrix isa UniformScaling || @boundscheck axes(mass_matrix) === axes(J) || _throwJMerror(J, mass_matrix)
  @inbounds if W_transform
    invdtgamma = inv(dtgamma)
    if MT <: UniformScaling
      λ = -mass_matrix.λ
      W = J + (λ * invdtgamma)*I
    else
      W = muladd(-mass_matrix, invdtgamma, J)
    end
  else
    if MT <: UniformScaling
      λ = -mass_matrix.λ
      W = dtgamma*J + λ*I
    else
      W = muladd(dtgamma, J, -mass_matrix)
    end
  end
  return W
end

function calc_W!(W, integrator, nlsolver::Union{Nothing,AbstractNLSolver}, cache, dtgamma, repeat_step, W_transform=false)
  @unpack t,dt,uprev,u,f,p = integrator
  lcache = nlsolver === nothing ? cache : nlsolver.cache
  @unpack J = lcache
  alg = unwrap_alg(integrator, true)
  mass_matrix = integrator.f.mass_matrix
  is_compos = integrator.alg isa CompositeAlgorithm
  isnewton = alg isa NewtonAlgorithm

  # handle Wfact
  if W_transform && DiffEqBase.has_Wfact_t(f)
    f.Wfact_t(W, u, p, dtgamma, t)
    is_compos && (integrator.eigen_est = opnorm(LowerTriangular(W), Inf) + inv(dtgamma)) # TODO: better estimate
    return nothing
  elseif !W_transform && DiffEqBase.has_Wfact(f)
    f.Wfact(W, u, p, dtgamma, t)
    if is_compos
      opn = opnorm(LowerTriangular(W), Inf)
      integrator.eigen_est = (opn + one(opn)) / dtgamma # TODO: better estimate
    end
    return nothing
  end

  # check if we need to update J or W
  #W_dt = nlsolver === nothing ? dt : nlsolver.cache.W_dt # TODO: RosW
  #new_jac = do_newJ(integrator, alg, cache, repeat_step)
  #new_W = do_newW(integrator, nlsolver, new_jac, W_dt)
  new_jac, new_W = do_newJW(integrator, alg, cache.nlsolver, repeat_step)

  (new_jac && isdefined(lcache, :J_t)) && (lcache.J_t = t)

  # calculate W
  if W isa WOperator
    isnewton || DiffEqBase.update_coefficients!(W,uprev,p,t) # we will call `update_coefficients!` in NLNewton
    W.transform = W_transform; set_gamma!(W, dtgamma)
  else # concrete W using jacobian from `calc_J!`
    islin, isode = islinearfunction(integrator)
    islin ? (J = isode ? f.f : f.f1.f) : ( new_jac && (calc_J!(J, integrator, lcache)) )
    new_W && jacobian2W!(W, mass_matrix, dtgamma, J, W_transform)
  end
  if isnewton
    set_new_W!(nlsolver, new_W) && set_W_dt!(nlsolver, dt)
  end
  new_W && (integrator.destats.nw += 1)
  return nothing
end

function calc_W(integrator, cache, dtgamma, repeat_step, W_transform=false)
  @unpack t,uprev,p,f = integrator
  mass_matrix = integrator.f.mass_matrix
  isarray = uprev isa AbstractArray
  # calculate W
  is_compos = integrator.alg isa CompositeAlgorithm
  islin, isode = islinearfunction(integrator)
  if islin
    J = isode ? f.f : f.f1.f # unwrap the Jacobian accordingly
    W = WOperator(mass_matrix, dtgamma, J, false; transform=W_transform)
  elseif DiffEqBase.has_jac(f)
    J = f.jac(uprev, p, t)
    if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
      J = DiffEqArrayOperator(J)
    end
    W = WOperator(mass_matrix, dtgamma, J, false; transform=W_transform)
    integrator.destats.nw += 1
  else
    integrator.destats.nw += 1
    J = calc_J(integrator, cache)
    W_full = W_transform ? -mass_matrix*inv(dtgamma) + J :
                           -mass_matrix + dtgamma*J
    W = W_full isa Number ? W_full : lu(W_full)
  end
  (W isa WOperator && unwrap_alg(integrator, true) isa NewtonAlgorithm) && (W = DiffEqBase.update_coefficients!(W,uprev,p,t)) # we will call `update_coefficients!` in NLNewton
  is_compos && (integrator.eigen_est = isarray ? opnorm(J, Inf) : abs(J))
  W
end

function calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, W_transform)
  calc_tderivative!(integrator, cache, dtd1, repeat_step)
  nlsolver = nothing
  # we need to skip calculating `W` when a step is repeated
  repeat_step || calc_W!(cache.W, integrator, nlsolver, cache, dtgamma, repeat_step, W_transform)
  return nothing
end

# update W matrix (only used in Newton method)
update_W!(integrator, cache, dt, repeat_step) =
  update_W!(cache.nlsolver, integrator, cache, dt, repeat_step)

function update_W!(nlsolver::AbstractNLSolver, integrator, cache::OrdinaryDiffEqMutableCache, dt, repeat_step)
  if isnewton(nlsolver)
    calc_W!(get_W(nlsolver), integrator, nlsolver, cache, dt, repeat_step, true)
  end
  nothing
end

function update_W!(nlsolver::AbstractNLSolver, integrator, cache, dt, repeat_step)
  if isnewton(nlsolver)
    nlsolver.cache.W = calc_W(integrator, nlsolver.cache, dt, repeat_step, true)
  end
  nothing
end

function build_J_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits,::Val{IIP}) where IIP
  islin, isode = islinearfunction(f, alg)
  if f.jac_prototype isa DiffEqBase.AbstractDiffEqLinearOperator
    W = WOperator(f, dt, IIP)
    J = W.J
  elseif IIP && f.jac_prototype !== nothing
    J = similar(f.jac_prototype)
    W = similar(J)
  elseif islin || (!IIP && DiffEqBase.has_jac(f))
    J = islin ? (isode ? f.f : f.f1.f) : f.jac(uprev, p, t) # unwrap the Jacobian accordingly
    if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
      J = DiffEqArrayOperator(J)
    end
    W = WOperator(f.mass_matrix, dt, J, IIP)
  else
    J = false .* _vec(u) .* _vec(u)'
    W = if IIP
      similar(J)
    else
      W = if u isa StaticArray
      lu(J)
      elseif u isa Number
        u
      else
        LU{LinearAlgebra.lutype(uEltypeNoUnits)}(Matrix{uEltypeNoUnits}(undef, 0, 0),
                                                     Vector{LinearAlgebra.BlasInt}(undef, 0),
                                                     zero(LinearAlgebra.BlasInt))
      end
    end # end W
  end
  return J, W
end

build_uf(alg::Union{DAEAlgorithm,OrdinaryDiffEqAlgorithm},nf,t,p,::Val{true}) =
  DiffEqDiffTools.UJacobianWrapper(nf,t,p)
build_uf(alg::Union{DAEAlgorithm,OrdinaryDiffEqAlgorithm},nf,t,p,::Val{false}) =
  DiffEqDiffTools.UDerivativeWrapper(nf,t,p)
