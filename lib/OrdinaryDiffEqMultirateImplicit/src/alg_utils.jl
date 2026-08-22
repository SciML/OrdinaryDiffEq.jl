alg_order(alg::MREIL) = alg.order
isfsal(::MREIL) = false

# The linearly implicit base method differentiates the *fast* component, so the
# shared `nlsolve_f`/`islinearfunction` machinery must resolve a SplitFunction to
# `f1` rather than to the whole right-hand side.
issplit(::MREIL) = true

# Inherited from the Rosenbrock algorithm union, but the linearly implicit Euler
# base method has no time-derivative term, so `t` is never differentiated through.
SciMLBase.forwarddiffs_model_time(::MREIL) = false

# The empty-pattern rejection has to run before the generic AD preparation:
# `prepare_user_sparsity` writes the mass-matrix diagonal into a sparse pattern in
# place, so by `alg_cache` time an empty prototype has already been given a diagonal —
# and a diagonal-only coloring of a non-diagonal `f1` then folds row sums into the
# diagonal of `J`, which is even further from the truth than the `J = 0` the empty
# pattern asks for. The `invoke` hands the algorithm on to the generic method, so
# `prepare_ADType`/`prepare_user_sparsity` still run exactly as without this method;
# `alg_cache` keeps the same check as a backstop for entry points that skip
# `prepare_alg` (a scalar state, or a user-supplied `AutoSparse`).
function DiffEqBase.prepare_alg(alg::MREIL, u0::AbstractArray, p, prob)
    _mreil_check_ad_sparsity(prob.f)
    return invoke(
        DiffEqBase.prepare_alg,
        Tuple{OrdinaryDiffEqAdaptiveImplicitAlgorithm, AbstractArray, Any, Any},
        alg, u0, p, prob
    )
end

alg_order(::MRIGARKIRK21a) = 2
alg_order(::MRIGARKESDIRK34a) = 3
alg_order(::MRIGARKESDIRK46a) = 4
isfsal(::MRIGARKImplicitAlg) = false

nlsolve_f(f, ::MRIGARKImplicitAlg) = f isa SplitFunction ? f.f2 : f
