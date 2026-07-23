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

Policy for whether the [`AutoDePSpecialize`](@ref) path de-specializes `p`.
`isbits` structs are packed into a `RespecializeParams.OpaqueParams` (byte copy)
and non-`isbits` structs into a `RespecializeParams.OpaqueRef` (by reference);
`RespecializeParams.pack_auto` / `opaque_container_type` pick the container and
the wrapper signature carries it in the `p` slot.

Excluded (kept as-is): `SciMLBase.NullParameters` (already a uniform singleton)
and already-packed containers (idempotency, so a sensitivity adjoint that
re-concretizes an opaque problem does not nest wrappers).

De-specializing makes `sol.prob.p` an opaque container, so **`AutoDePSpecialize`
is a forward-solve latency tool, not for parameter-estimation workflows**:
reverse-mode AD of `p` (SciMLSensitivity) and symbolic parameter indexing
(`getp`/`setp`/observed) both read the concrete parameter object and are not
supported under it (reverse-mode support is a planned follow-up). Recover the
payload explicitly with `RespecializeParams.unpack(sol.prob.p, typeof(p))`.
"""
@inline should_opaque_p(p) = !(p isa SciMLBase.NullParameters) &&
    !(p isa RespecializeParams.OpaqueParams) && !(p isa RespecializeParams.OpaqueRef)

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
