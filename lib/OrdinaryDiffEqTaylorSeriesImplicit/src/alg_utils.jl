alg_extrapolates(::ImplicitTaylor) = true

# not used. just to make API satisfied
alg_autodiff(::ImplicitTaylor) = AutoForwardDiff()

function alg_order(alg::ImplicitTaylor{P}) where {P}
    μ = alg.μ
    μ̄ = 1 - μ
    P_plus_1_coeff = if isodd(P)
        μ̄^(P + 1) - μ^(P + 1)
    else
        μ̄^(P + 1) + μ^(P + 1)
    end
    if P_plus_1_coeff ≈ 0
        # when the function is real, we can get an additional order by projection
        return alg.real_function ? P + 2 : P + 1
    else
        return P
    end
end

alg_adaptive_order(::ImplicitTaylor) = 0
