"""
    JVPCache{T} <: AbstractSciMLOperator{T}

JVPCache provides a JVP operator wrapper for performing the DifferentiationInterface pushforward operation.

### Constructor

```julia
JVPCache(f::SciMLBase.AbstractDiffEqFunction, du, u, p, t; autodiff)
```

JVPCache construction builds a DifferentiationInterface "prep" object using `prepare_pushforward!`. The "prep" object is used
when applying the operator.

### Computing the JVP

Computing the JVP is done with the DifferentiationInterface function `pushforward!`, which takes advantage of the preparation done upon construction.

### Counting

`njvps` tallies the products applied since the count was last drained by
`drain_jvp_count!`, which turns them into the `stats.nf` they cost. The tally
lives here, on the operator that performs the products, because that is the only place
that sees all of them: a Krylov solve applies the operator more times than its reported
iteration count (a warm start and the initial residual each cost one), and one operator
can be shared by several `W`s.
"""
@concrete mutable struct JVPCache{T} <: SciMLOperators.AbstractSciMLOperator{T}
    jvp_op::Any
    f::Any
    du::Any
    u::Any
    p::Any
    t::Any
    njvps::Int
end

SciMLBase.isinplace(::JVPCache) = true
ArrayInterface.can_setindex(::JVPCache) = false
# A JVPCache is genuinely matrix-free: it only implements `mul!` (the DI pushforward), with
# no `convert(AbstractMatrix, ·)`. Report that honestly so a `WOperator` wrapping one is
# itself matrix-free and consumers (e.g. NonlinearSolve) route it to a matrix-free linear
# solver instead of trying to factorize it.
SciMLOperators.isconvertible(::JVPCache) = false
function ArrayInterface.restructure(y::JVPCache, x::JVPCache)
    @assert size(y) == size(x) "cannot restructure operators. ensure their sizes match."
    return x
end

function ConstructionBase.constructorof(::Type{<:JVPCache{T}}) where {T}
    return JVPCache{T}
end

Base.size(J::JVPCache) = (length(J.u), length(J.u))

function JVPCache(f::SciMLBase.AbstractDiffEqFunction, du, u, p, t; autodiff)
    jvp_op = prepare_jvp(f, du, u, p, t, autodiff)
    return JVPCache{eltype(du)}(jvp_op, f, du, u, p, t, 0)
end

function (op::JVPCache)(Jv, v, u, p, t)
    op.jvp_op(Jv, v, u, p, t)
    op.njvps += 1
    return Jv
end

function LinearAlgebra.mul!(
        Jv::AbstractArray, J::JVPCache, v::AbstractArray
    )
    J.jvp_op(Jv, v, J.u, J.p, J.t)
    J.njvps += 1
    return Jv
end

"""
    rhs_evals_per_jvp(ad) -> Int

Right-hand-side evaluations one matrix-free Jacobian-vector product costs under autodiff
backend `ad`.

`stats.nf` counts calls to `f`, not cost-equivalent work. A finite-difference JVP calls
`f` twice — once at the linearization point and once at the perturbed point — while a
dual-number JVP calls it once, on a dual input that costs more per call but is still one
call. Reporting the finite-difference product as one evaluation to make the two look
comparable would put a cost model into a counter that nothing else in the solver treats
as one.

The finite-difference count drops to one if the base evaluation is ever cached across the
products of a single linear solve (`FiniteDiff.finite_difference_jvp!` takes an `f_in`
argument for exactly that, which DifferentiationInterface's `pushforward!` does not
currently supply).
"""
function rhs_evals_per_jvp(ad)
    ad = ad isa AutoSparse ? ADTypes.dense_ad(ad) : ad
    return (ad isa AutoFiniteDiff || ad isa ADTypes.AutoFiniteDifferences) ? 2 : 1
end

"""
    drain_jvp_count!(integrator, alg, W) -> nothing

Add the Jacobian-vector products accumulated in `W`'s JVP operator to `integrator.stats.nf`
and reset the tally, so each product is counted once no matter how many `W`s share the
operator.

Nothing is counted when no `JVPCache` is involved: a `W` whose Jacobian is a concrete
matrix or a user-supplied `MatrixOperator` applies a stored matrix and evaluates `f` zero
times.
"""
function drain_jvp_count!(integrator, alg, W)
    J = jvp_counter(W)
    (J === nothing || !(integrator isa SciMLBase.DEIntegrator)) && return nothing
    n = J.njvps
    J.njvps = 0
    n == 0 && return nothing
    OrdinaryDiffEqCore.increment_nf!(
        integrator.stats, rhs_evals_per_jvp(alg_autodiff(alg)) * n
    )
    return nothing
end

"""
    jvp_counter(W) -> JVPCache or nothing

The `JVPCache` whose products `W` applies, or `nothing` when `W` costs no RHS evaluations
per product. A `WOperator` uses its `jacvec` when it has one and its `J` otherwise, which
is what `LinearAlgebra.mul!` does with it.
"""
jvp_counter(J::JVPCache) = J
jvp_counter(W::WOperator) = jvp_counter(W.jacvec === nothing ? W.J : W.jacvec)
jvp_counter(::Any) = nothing

# helper functions

function prepare_jvp(f::SciMLBase.AbstractDiffEqFunction, du, u, p, t, autodiff)
    SciMLBase.has_jvp(f) && return f.jvp
    autodiff = autodiff isa AutoSparse ? ADTypes.dense_ad(autodiff) : autodiff
    @assert DI.check_inplace(autodiff) "AD backend $(autodiff) doesn't support in-place problems."
    di_prep = DI.prepare_pushforward(
        f, du, autodiff, u, (u,), DI.ConstantOrCache(p), DI.Constant(t)
    )
    return (
        Jv,
        v,
        u,
        p,
        t,
    ) -> DI.pushforward!(
        f, du, (reshape(Jv, size(du)),), di_prep, autodiff, u,
        (reshape(v, size(u)),), DI.ConstantOrCache(p), DI.Constant(t)
    )
end

function SciMLOperators.update_coefficients!(J::JVPCache, u, p, t)
    J.u = u
    J.p = p
    return J.t = t
end

function resize_JVPCache!(J::JVPCache, f, du, u, p, t, autodiff)
    J.jvp_op = prepare_jvp(f, du, u, p, t, autodiff)
    J.du = du
    return update_coefficients!(J, u, p, t)
end

function resize_JVPCache!(J::JVPCache, f, du, u, autodiff)
    J.jvp_op = prepare_jvp(f, du, u, J.p, J.t, autodiff)
    J.du = du
    return update_coefficients!(J, u, J.p, J.t)
end
