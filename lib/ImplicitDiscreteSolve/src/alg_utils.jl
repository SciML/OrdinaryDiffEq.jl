function SciMLBase.isautodifferentiable(alg::IDSolve)
    return true
end
function SciMLBase.allows_arbitrary_number_types(alg::IDSolve)
    return true
end
function SciMLBase.allowscomplex(alg::IDSolve)
    return true
end

SciMLBase.isdiscrete(alg::IDSolve) = true

isfsal(alg::IDSolve) = false
alg_order(alg::IDSolve) = 1
beta2_default(alg::IDSolve) = 0
beta1_default(alg::IDSolve, beta2) = 0

dt_required(alg::IDSolve) = false
isdiscretealg(alg::IDSolve) = true

isadaptive(alg::IDSolve) = true

# Preserve `u === nothing` through __init so alg_cache(alg, ::Nothing, ...) fires
# and we skip building a NonlinearProblem over an empty state.
allows_null_u0(alg::IDSolve) = true
