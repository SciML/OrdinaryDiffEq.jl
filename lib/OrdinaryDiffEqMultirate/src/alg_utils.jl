alg_order(alg::MREEF) = alg.order
isfsal(::MREEF) = false

function prepare_alg(alg::MREEF, u0::AbstractArray, p, prob)
    return alg
end

alg_order(alg::MRAB) = alg.k
isfsal(::MRAB) = false

function prepare_alg(alg::MRAB, u0::AbstractArray, p, prob)
    1 <= alg.k <= 5 || throw(ArgumentError("MRAB: `k` must satisfy 1 ≤ k ≤ 5"))
    alg.m >= 1 || throw(ArgumentError("MRAB: `m` must be ≥ 1"))
    return alg
end

alg_order(::MRIGARKERK22a) = 2
isfsal(::MRIGARKERK22a) = false

function prepare_alg(alg::MRIGARKERK22a, u0::AbstractArray, p, prob)
    alg.m >= 1 || throw(ArgumentError("MRIGARKERK22a: `m` must be ≥ 1"))
    return alg
end
