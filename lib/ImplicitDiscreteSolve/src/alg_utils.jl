function SciMLBase.isautodifferentiable(alg::IDSolve)
    true
end
function SciMLBase.allows_arbitrary_number_types(alg::IDSolve)
    true
end
function SciMLBase.allowscomplex(alg::IDSolve)
    true
end

SciMLBase.isdiscrete(alg::IDSolve) = true

isfsal(alg::IDSolve) = false
alg_order(alg::IDSolve) = 0
beta2_default(alg::IDSolve) = 0
beta1_default(alg::IDSolve, beta2) = 0

dt_required(alg::IDSolve) = false
isdiscretealg(alg::IDSolve) = true
