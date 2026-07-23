# AutoDePSpecialize ODE-path wrapping.
#
# The pieces this path uses that are shared, not defined here:
#   * `AutoDePSpecialize` marker — SciMLBase.
#   * `OpaqueParams`/`OpaqueRef`/`OpaqueVoid`, `pack_auto`, and the
#     `wrap_void_opaque` installer — RespecializeParams. `promote_f` calls
#     `RespecializeParams.wrap_void_opaque(ff, P, (natural_sig,))` to build the
#     de-specialized `FunctionWrappersWrapper` for the RHS / tgrad / jac.
#
# `should_opaque_p` (the applicability policy) is defined here rather than in
# SciMLBase so this does not depend on an unreleased SciMLBase symbol. It is a
# trivial predicate; a solver stack that wants to share it (e.g. NonlinearSolve)
# can promote it to SciMLBase in a follow-up.
#
# The solution returned from `solve` keeps `prob.p` as the packed
# `OpaqueParams`/`OpaqueRef` container — it must stay that way because dense/lazy
# interpolation re-invokes the wrapped RHS with `sol.prob.p`. Recover the
# original payload with `RespecializeParams.unpack(sol.prob.p, P)`.

"""
    should_opaque_p(p) :: Bool

Policy for whether the [`AutoDePSpecialize`](@ref) path de-specializes `p`: true
for an `isbits` non-`NullParameters` payload, false otherwise. `NullParameters`
is already a uniform singleton, and non-`isbits` payloads fall back to the plain
`AutoSpecialize` wrapping.
"""
@inline should_opaque_p(p) = isbits(p) && !(p isa SciMLBase.NullParameters)

"""
    wrapfun_iip_opaque(ff, ::Type{P}, inputs, ::Val{CS}) -> FunctionWrappersWrapper

Opaque-aware analogue of [`wrapfun_iip`](@ref): builds the four-variant
(plain / Jacobian-`Dual` / time-`Dual` / both) wrapper with the opaque container
substituted in the `p` slot. `inputs` is the original `(u0, u0, p, t)` tuple.
Defined in `DiffEqBaseForwardDiffExt`; only called from the `Val{true}`
(ForwardDiff-using) `promote_f` path, which requires the extension anyway. The
non-ForwardDiff path uses `RespecializeParams.wrap_void_opaque` directly so that
the wrapper type does not depend on whether ForwardDiff is loaded.
"""
function wrapfun_iip_opaque end
