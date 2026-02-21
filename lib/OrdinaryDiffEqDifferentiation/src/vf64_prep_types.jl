# Concrete DI prep type aliases for the VF64 FiniteDiff path.
#
# FiniteDiff jacobian preps are function-independent (unlike ForwardDiff which encodes
# the function type in its Tag). This means we can compute the concrete type at module
# load time and hardcode it in VF64 structs, eliminating JCType as a type parameter.
#
# Note: FiniteDiff derivative preps (grad_config) DO encode the function type, so
# GCType cannot be eliminated for Rosenbrock caches.
# ForwardDiff preps encode both function Tag and chunk size, so JCType cannot be
# eliminated for the ForwardDiff path either.

let
    n = 2
    _du = zeros(n)
    _u = zeros(n)
    _uf_dummy = (y, x) -> nothing

    _fd_fwd = DI.prepare_jacobian(
        _uf_dummy, _du, AutoFiniteDiff(; dir = 1), _u, strict = Val(false))
    _fd_rev = DI.prepare_jacobian(
        _uf_dummy, _du, AutoFiniteDiff(; dir = -1), _u, strict = Val(false))
    global const _JacPrepFiniteDiffForward = typeof(_fd_fwd)
    global const _JacPrepFiniteDiffReverse = typeof(_fd_rev)
    global const _JacConfigFiniteDiff = Tuple{_JacPrepFiniteDiffForward, _JacPrepFiniteDiffReverse}
end
