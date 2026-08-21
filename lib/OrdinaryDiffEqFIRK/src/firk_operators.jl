"""
    ComplexWOperator{T, MMType, GType, JType, JVType, VType} <: AbstractSciMLOperator{T}

Lazy representation of the complex-shifted FIRK stage matrix

```math
W = -\\frac{MM}{\\gamma} + J
```

with complex `gamma` but a *real* mass matrix `MM` and a *real* Jacobian `J`. It is the
complex counterpart of `SciMLOperators.WOperator` and follows the same `gamma` convention,
so `-γ_c MM + J` is obtained with `gamma = inv(γ_c)`.

Because `J` is real, the complex product is assembled from two real Jacobian-vector
products, `J (a + bi) = Ja + i Jb`. That is what allows a matrix-free `J` — a real
Jacobian-vector-product operator, which cannot be applied to a complex vector — to back
the complex-conjugate stage pairs of the FIRK W-transform.

# Arguments

  - `mass_matrix`: `UniformScaling` or real matrix `MM`.
  - `gamma`: complex scalar; the prototype value fixes the operator's element type.
  - `J`: real Jacobian, either a matrix or an `AbstractSciMLOperator`.
  - `u`: prototype state, used to allocate the real work buffers.
  - `jacvec`: optional real Jacobian-vector-product operator applied in place of `J`.
"""
mutable struct ComplexWOperator{T, MMType, GType, JType, JVType, VType} <:
    AbstractSciMLOperator{T}
    mass_matrix::MMType
    gamma::GType
    J::JType
    jacvec::JVType
    _rev::VType
    _imv::VType
    _reJv::VType
    _imJv::VType
end

function ComplexWOperator(mass_matrix, gamma::Number, J, u, jacvec = nothing)
    v = _vec(u)
    T = complex(promote_type(eltype(v), real(typeof(gamma))))
    return ComplexWOperator{
        T, typeof(mass_matrix), T, typeof(J), typeof(jacvec), typeof(v),
    }(
        mass_matrix, convert(T, gamma), J, jacvec,
        zero(v), zero(v), zero(v), zero(v)
    )
end

Base.eltype(::ComplexWOperator{T}) where {T} = T
Base.size(W::ComplexWOperator) = size(W.J)
Base.size(W::ComplexWOperator, d::Integer) = d <= 2 ? size(W)[d] : 1
SciMLBase.isinplace(::ComplexWOperator) = true
SciMLOperators.isconvertible(::ComplexWOperator) = false

_jacobian_operator(W::ComplexWOperator) = W.jacvec === nothing ? W.J : W.jacvec

function SciMLOperators.update_coefficients!(
        W::ComplexWOperator, u = nothing, p = nothing, t = nothing; gamma = nothing
    )
    if u !== nothing && p !== nothing && t !== nothing
        SciMLOperators.update_coefficients!(W.J, u, p, t)
        W.jacvec === nothing || SciMLOperators.update_coefficients!(W.jacvec, u, p, t)
    end
    gamma === nothing || (W.gamma = gamma)
    return W
end

function LinearAlgebra.mul!(Y::AbstractVector, W::ComplexWOperator, B::AbstractVector)
    (; _rev, _imv, _reJv, _imJv) = W
    Jop = _jacobian_operator(W)
    @. _rev = real(B)
    @. _imv = imag(B)
    mul!(_reJv, Jop, _rev)
    mul!(_imJv, Jop, _imv)
    @. Y = complex(_reJv, _imJv)
    a = -inv(W.gamma)
    if W.mass_matrix isa UniformScaling
        @. Y += (a * W.mass_matrix.λ) * B
    else
        mul!(_reJv, W.mass_matrix, _rev)
        mul!(_imJv, W.mass_matrix, _imv)
        @. Y += a * complex(_reJv, _imJv)
    end
    return Y
end

function LinearAlgebra.mul!(Y::AbstractArray, W::ComplexWOperator, B::AbstractArray)
    mul!(vec(Y), W, vec(B))
    return Y
end

function Base.:*(W::ComplexWOperator, x::AbstractVector)
    return mul!(similar(x, promote_type(eltype(W), eltype(x))), W, x)
end

const LazyW = Union{SciMLOperators.WOperator, ComplexWOperator}

"""
    is_lazy_W(W) -> Bool

Whether the FIRK stage matrix `W` is applied lazily rather than stored as a concrete
matrix assembled entrywise from `J`.
"""
is_lazy_W(W) = W isa LazyW

"""
    set_W_gamma!(W, γ)

Point the lazily-applied FIRK stage matrix `W` at `-γ MM + J` for the current step. Both
`SciMLOperators.WOperator` (real stage) and [`ComplexWOperator`](@ref) store the reciprocal
shift, hence the `inv`.

This is O(1) and no factorization follows it, so the callers re-form a lazy `W` on every
step rather than reusing a stale shift the way the concrete path does — there is nothing
to amortize, and a current `W` is a better Newton matrix.
"""
set_W_gamma!(W::LazyW, γ) = (W.gamma = inv(γ); W)

"""
    firk_real_W(alg, J, W) -> W

Stage matrix for the real FIRK stage. `build_J_W` returns a `WOperator` whenever a
Jacobian-vector product is available, which includes `concrete_jac = true` with a
factorization `linsolve`. FIRK assembles its own shifted matrices entrywise, so a lazy `W`
over a concrete `J` is only useful when the linear solver either is itself matrix-free or
consumes the split form directly (`LHLFactorization`, which reduces `J` once and absorbs
each stage shift in O(n²)); otherwise fall back to a concrete matrix.
"""
function firk_real_W(alg, J, W)
    (is_lazy_W(W) && !(J isa AbstractSciMLOperator)) || return W
    alg.linsolve isa LinearSolve.LHLFactorization && return W
    (alg.linsolve === nothing || LinearSolve.needs_concrete_A(alg.linsolve)) || return W
    return recursivefill!(similar(J), false)
end

"""
    build_complex_W(alg, f, u, J, W1) -> W

Stage matrix for a complex-conjugate FIRK stage pair: a concrete complex matrix when the
real stage matrix `W1` is concrete, and a lazily-applied [`ComplexWOperator`](@ref) when
`W1` is itself an operator (matrix-free Jacobian and/or Krylov `linsolve`).
"""
function build_complex_W(alg, f, u, J, W1)
    is_lazy_W(W1) || return recursivefill!(similar(J, Complex{eltype(W1)}), false)
    if alg.linsolve isa LinearSolve.LHLFactorization && !(J isa AbstractSciMLOperator)
        # Same split `-MM/gamma + J` shape as the real stage, only with a complex `gamma`
        # slot. The Jacobian object is shared with the real stage — one `calc_J!` serves
        # every stage matrix — while inside the linear solver the reduction stays real and
        # only the shift goes complex.
        return SciMLOperators.WOperator{true}(
            f.mass_matrix, complex(one(eltype(W1))), J, complex.(_vec(u))
        )
    end
    if alg.linsolve === nothing || LinearSolve.needs_concrete_A(alg.linsolve)
        error(
            "$(nameof(typeof(alg))) got a non-concrete (operator) Jacobian, which can only " *
                "be used with a matrix-free linear solver, but `linsolve = $(alg.linsolve)` " *
                "requires a concrete matrix. Either pass a Krylov solver (e.g. " *
                "`linsolve = KrylovJL_GMRES()`), or give the problem a concrete " *
                "`jac_prototype` so that `W` can be assembled and factorized."
        )
    end
    jacvec = W1 isa SciMLOperators.WOperator ? W1.jacvec : nothing
    return ComplexWOperator(f.mass_matrix, complex(one(eltype(W1))), J, u, jacvec)
end


"""
    firk_new_J!(J, W, integrator, cache)

Refresh the Jacobian for a new step: recompute a concrete `J` in place, or move a
matrix-free Jacobian's linearization point to `(uprev, p, t)`. `W` carries the
Jacobian-vector-product operator when one exists alongside a concrete `J`.
"""
function firk_new_J!(J, W, integrator, cache, extra_Ws...)
    (; uprev, p, t) = integrator
    if J isa AbstractSciMLOperator
        SciMLOperators.update_coefficients!(J, uprev, p, t)
    else
        calc_J!(J, integrator, cache)
    end
    jacvec = W isa LazyW ? W.jacvec : nothing
    if jacvec !== nothing && jacvec !== J
        SciMLOperators.update_coefficients!(jacvec, uprev, p, t)
    end
    # Every stage matrix viewing this Jacobian must hear that it moved; a linear solver
    # caching a factorization of `J` across steps (LHLFactorization) has no other way to
    # tell an in-place refresh from the unchanged object.
    _mark_J_moved!(W)
    foreach(_mark_J_moved!, extra_Ws)
    return nothing
end

_mark_J_moved!(W::SciMLOperators.WOperator) = SciMLOperators.mark_jacobian_updated!(W)
_mark_J_moved!(Ws::AbstractVector) = foreach(_mark_J_moved!, Ws)
_mark_J_moved!(::Any) = nothing

# A `ComplexWOperator` applies the same real JVP operator its `W` does, so the products it
# performs are already on that operator's tally and `drain_jvp_count!` needs nothing else.
jvp_counter(W::ComplexWOperator) = jvp_counter(_jacobian_operator(W))
