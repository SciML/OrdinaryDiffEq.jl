function SciMLBase.isautodifferentiable(alg::FunctionMap)
    true
end
function SciMLBase.allows_arbitrary_number_types(alg::FunctionMap)
    true
end
function SciMLBase.allowscomplex(alg::FunctionMap)
    true
end

SciMLBase.isdiscrete(alg::FunctionMap) = true

isfsal(alg::FunctionMap) = false

alg_order(alg::FunctionMap) = 0

beta2_default(alg::FunctionMap) = 0

beta1_default(alg::FunctionMap, beta2) = 0

function FunctionMap_scale_by_time(alg::FunctionMap{scale_by_time}) where {scale_by_time}
    scale_by_time
end