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

alg_order(::Union{MRIGARKERK22a, MRIGARKERK22b}) = 2
alg_order(::MRIGARKERK33a) = 3
isfsal(::Union{MRIGARKERK22a, MRIGARKERK22b, MRIGARKERK33a}) = false

function prepare_alg(
        alg::Union{MRIGARKERK22a, MRIGARKERK22b, MRIGARKERK33a}, u0::AbstractArray, p, prob
    )
    alg.m >= 1 || throw(ArgumentError("$(nameof(typeof(alg))): `m` must be ≥ 1"))
    return alg
end

alg_order(::MIS) = 2
isfsal(::MIS) = false

function prepare_alg(alg::MIS, u0::AbstractArray, p, prob)
    alg.m >= 1 || throw(ArgumentError("MIS: `m` must be ≥ 1"))
    return alg
end
