alg_order(alg::MREIL) = alg.order
isfsal(::MREIL) = false

# The linearly implicit base method differentiates the *fast* component, so the
# shared `nlsolve_f`/`islinearfunction` machinery must resolve a SplitFunction to
# `f1` rather than to the whole right-hand side.
issplit(::MREIL) = true

# Inherited from the Rosenbrock algorithm union, but the linearly implicit Euler
# base method has no time-derivative term, so `t` is never differentiated through.
SciMLBase.forwarddiffs_model_time(::MREIL) = false

alg_order(::MRIGARKIRK21a) = 2
alg_order(::MRIGARKESDIRK34a) = 3
alg_order(::MRIGARKESDIRK46a) = 4
isfsal(::MRIGARKImplicitAlg) = false

function prepare_alg(alg::MRIGARKImplicitAlg, u0::AbstractArray, p, prob)
    alg.m >= 1 || throw(ArgumentError("$(nameof(typeof(alg))): `m` must be ≥ 1"))
    return alg
end

nlsolve_f(f, ::MRIGARKImplicitAlg) = f isa SplitFunction ? f.f2 : f
