function SciMLBase.isautodifferentiable(alg::FunctionMap)
    return true
end
function SciMLBase.allows_arbitrary_number_types(alg::FunctionMap)
    return true
end
function SciMLBase.allowscomplex(alg::FunctionMap)
    return true
end

SciMLBase.isdiscrete(alg::FunctionMap) = true

isfsal(alg::FunctionMap) = false

alg_order(alg::FunctionMap) = 0

beta2_default(alg::FunctionMap) = 0

beta1_default(alg::FunctionMap, beta2) = 0

function FunctionMap_scale_by_time(alg::FunctionMap{scale_by_time}) where {scale_by_time}
    return scale_by_time
end

dt_required(alg::FunctionMap) = false
isdiscretealg(alg::FunctionMap) = true
