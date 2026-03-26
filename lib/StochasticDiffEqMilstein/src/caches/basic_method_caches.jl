struct RKMilGeneralConstantCache{JalgType, dWType, dZType} <: StochasticDiffEqConstantCache
    Jalg::JalgType
    # For rejection handling: store original dW, dZ, dt from the first attempt
    # so sub-interval iterated integrals can be computed from the same Fourier coefficients
    _dW_orig::Base.RefValue{dWType}
    _dZ_orig::Base.RefValue{dZType}
    _dt_orig::Base.RefValue{Float64}
end

@cache struct RKMilGeneralCache{uType, rateType, rateNoiseType, JalgType, dWType, dZType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    du1::rateType
    du2::rateType
    K::uType
    L::rateNoiseType
    mil_correction::uType
    ggprime::rateNoiseType
    Jalg::JalgType
    _dW_orig::Base.RefValue{dWType}
    _dZ_orig::Base.RefValue{dZType}
    _dt_orig::Base.RefValue{Float64}
end

function alg_cache(
        alg::RKMilGeneral, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    Jalg = get_Jalg(DeltaW, dt, prob, alg)
    dW_ref = Ref(copy(DeltaW))
    dZ_ref = DeltaZ === nothing ? Ref{Nothing}(nothing) : Ref(copy(DeltaZ))
    return RKMilGeneralConstantCache(Jalg, dW_ref, dZ_ref, Ref(0.0))
end

function alg_cache(
        alg::RKMilGeneral, prob, u, DeltaW, DeltaZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    K = zero(u)
    L = zero(noise_rate_prototype)
    mil_correction = zero(u)
    ggprime = zero(noise_rate_prototype)
    Jalg = get_Jalg(DeltaW, dt, prob, alg)
    dW_ref = Ref(copy(DeltaW))
    dZ_ref = DeltaZ === nothing ? Ref{Nothing}(nothing) : Ref(copy(DeltaZ))
    return RKMilGeneralCache(
        u, uprev, tmp, du1, du2, K, L, mil_correction, ggprime, Jalg,
        dW_ref, dZ_ref, Ref(0.0)
    )
end
