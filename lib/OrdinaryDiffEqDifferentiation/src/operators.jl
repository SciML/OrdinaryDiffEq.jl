"""
    JVPCache{T} <: AbstractSciMLOperator{T}

JVPCache provides a JVP operator wrapper for performing the DifferentiationInterface pushforward operation.

### Constructor

```julia
JVPCache(f::SciMLBase.AbstractDiffEqFunction, du, u, p, t; autodiff)
```

JVPCache construction builds a DifferentiationInterface "prep" object using `prepare_pushforward`. The "prep" object is used
when applying the operator.

### Computing the JVP

Computing the JVP is done with the DifferentiationInterface function `pushforward!`, which takes advantage of the preparation done upon construction.

The linearization point only moves on `update_coefficients!`, so the prep is additionally
specialized to it via `prepare_pushforward_same_point` (see `sync_jvp_point!`). A
finite-difference backend then evaluates `f` at the point once and reuses it across every
product taken there, instead of re-evaluating it per product.

### Counting

`njvps` tallies the products applied and `nbase_evals` the RHS evaluations spent fixing the
linearization point, both since the count was last drained by `drain_jvp_count!`, which
turns them into the `stats.nf` they cost. The tally lives here, on the operator that
performs the products, because that is the only place that sees all of them: a Krylov solve
applies the operator more times than its reported iteration count (a warm start and the
initial residual each cost one), and one operator can be shared by several `W`s.
"""
@concrete mutable struct JVPCache{T} <: SciMLOperators.AbstractSciMLOperator{T}
    jvp_op::Any
    f::Any
    du::Any
    u::Any
    p::Any
    t::Any
    njvps::Int
    nbase_evals::Int
    stale_point::Bool
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
    return JVPCache{eltype(du)}(jvp_op, f, du, u, p, t, 0, 0, true)
end

function (op::JVPCache)(Jv, v, u, p, t)
    update_coefficients!(op, u, p, t)
    return LinearAlgebra.mul!(Jv, op, v)
end

function LinearAlgebra.mul!(
        Jv::AbstractArray, J::JVPCache, v::AbstractArray
    )
    sync_jvp_point!(J)
    J.jvp_op(Jv, v, J.u, J.p, J.t)
    J.njvps += 1
    return Jv
end

"""
    sync_jvp_point!(J::JVPCache) -> nothing

Fix the operator's linearization point at `J.u`, `J.p`, `J.t`, so that a backend able to
reuse the base evaluation `f(u)` across products spends it once per point rather than once
per product.

Deferred to the first product after `update_coefficients!` rather than done there, so that a
point which is set but never applied costs nothing, and so that the two
`update_coefficients!` calls a `WOperator` makes when its `J` and its `jacvec` are the same
operator still fix the point once.
"""
function sync_jvp_point!(J::JVPCache)
    J.stale_point || return nothing
    J.jvp_op, n = set_jvp_point!(J.jvp_op, J.u, J.p, J.t)
    J.stale_point = false
    J.nbase_evals += n
    return nothing
end

"""
    rhs_evals_per_jvp(ad) -> Int

Right-hand-side evaluations one matrix-free Jacobian-vector product costs under autodiff
backend `ad`.

`stats.nf` counts calls to `f`, not cost-equivalent work. A dual-number JVP calls `f` once,
on a dual input that costs more per call but is still one call. A finite-difference JVP
evaluates only the perturbed point, its base evaluation having been paid once for the whole
linearization point and counted by `base_evals_per_jvp_point`. `AutoFiniteDifferences` is
the one backend `DifferentiationInterface` gives no such reuse, so it keeps its second
evaluation.
"""
function rhs_evals_per_jvp(ad)
    ad = ad isa AutoSparse ? ADTypes.dense_ad(ad) : ad
    return ad isa ADTypes.AutoFiniteDifferences ? 2 : 1
end

"""
    base_evals_per_jvp_point(ad) -> Int

RHS evaluations `set_jvp_point!` costs under autodiff backend `ad`, paid once per
linearization point and amortized over every product taken at it.

`AutoFiniteDiff` is the only backend whose
`DifferentiationInterface.prepare_pushforward_same_point` evaluates `f`; for the rest it is
a no-op.
"""
function base_evals_per_jvp_point(ad)
    ad = ad isa AutoSparse ? ADTypes.dense_ad(ad) : ad
    return ad isa AutoFiniteDiff ? 1 : 0
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
    n, nbase = J.njvps, J.nbase_evals
    J.njvps, J.nbase_evals = 0, 0
    (n == 0 && nbase == 0) && return nothing
    OrdinaryDiffEqCore.increment_nf!(
        integrator.stats, rhs_evals_per_jvp(alg_autodiff(alg)) * n + nbase
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

"""
    DIPushforward(f, du, backend, prep)

Callable `(Jv, v, u, p, t)` applying `DifferentiationInterface.pushforward!` with a prep
object, which `set_jvp_point!` respecializes to each linearization point. Held immutably
because a backend may answer `prepare_pushforward_same_point` with a prep of a different
type than the one it was handed.
"""
struct DIPushforward{F, DU, B, P}
    f::F
    du::DU
    backend::B
    prep::P
end

function (op::DIPushforward)(Jv, v, u, p, t)
    DI.pushforward!(
        op.f, op.du, (reshape(Jv, size(op.du)),), op.prep, op.backend, u,
        (reshape(v, size(u)),), DI.ConstantOrCache(p), DI.Constant(t)
    )
    return Jv
end

"""
    set_jvp_point!(jvp_op, u, p, t) -> (jvp_op, nf)

Fix `jvp_op`'s linearization point at `(u, p, t)`, returning the operator to apply from now
on and the RHS evaluations doing so cost.

A user-supplied `f.jvp` takes its point per product and has nothing to fix.
"""
set_jvp_point!(jvp_op, u, p, t) = (jvp_op, 0)

function set_jvp_point!(op::DIPushforward, u, p, t)
    prep = DI.prepare_pushforward_same_point(
        op.f, op.du, op.prep, op.backend, u, (u,),
        DI.ConstantOrCache(p), DI.Constant(t)
    )
    # `AutoFiniteDiff` reuses the prep it was handed, so the common path stays allocation-free.
    new_op = prep === op.prep ? op : DIPushforward(op.f, op.du, op.backend, prep)
    return new_op, base_evals_per_jvp_point(op.backend)
end

function prepare_jvp(f::SciMLBase.AbstractDiffEqFunction, du, u, p, t, autodiff)
    SciMLBase.has_jvp(f) && return f.jvp
    autodiff = autodiff isa AutoSparse ? ADTypes.dense_ad(autodiff) : autodiff
    @assert DI.check_inplace(autodiff) "AD backend $(autodiff) doesn't support in-place problems."
    di_prep = DI.prepare_pushforward(
        f, du, autodiff, u, (u,), DI.ConstantOrCache(p), DI.Constant(t)
    )
    return DIPushforward(f, du, autodiff, di_prep)
end

function SciMLOperators.update_coefficients!(J::JVPCache, u, p, t)
    J.u = u
    J.p = p
    J.t = t
    J.stale_point = true
    return J.t
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
