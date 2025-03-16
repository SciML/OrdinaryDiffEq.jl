function SciMLBase.isautodifferentiable(alg::SimpleIDSolve)
    true
end
function SciMLBase.allows_arbitrary_number_types(alg::SimpleIDSolve)
    true
end
function SciMLBase.allowscomplex(alg::SimpleIDSolve)
    true
end

SciMLBase.isdiscrete(alg::SimpleIDSolve) = true

isfsal(alg::SimpleIDSolve) = false
alg_order(alg::SimpleIDSolve) = 0
beta2_default(alg::SimpleIDSolve) = 0
beta1_default(alg::SimpleIDSolve, beta2) = 0

dt_required(alg::SimpleIDSolve) = false
isdiscretealg(alg::SimpleIDSolve) = true
