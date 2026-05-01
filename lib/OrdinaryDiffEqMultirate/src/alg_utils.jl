alg_order(alg::MREEF) = alg.order
isfsal(::MREEF) = false

function prepare_alg(alg::MREEF, u0::AbstractArray, p, prob)
    return alg
end
