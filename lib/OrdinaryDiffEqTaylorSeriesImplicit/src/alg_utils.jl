alg_extrapolates(::ImplicitTaylor) = true

# not used. just to make API satisfied
alg_autodiff(::ImplicitTaylor) = AutoForwardDiff()

is_mu_taylor(::ImplicitTaylor{P, Q}) where {P, Q} = Q == 0

is_taylor_pade(::ImplicitTaylor{P, Q}) where {P, Q} = Q != 0

function alg_order(alg::ImplicitTaylor{P, Q}) where {P, Q}
    if is_mu_taylor(alg) # μ-Taylor method
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
    else # Taylor-Padé method
        return P + Q
    end
end

function alg_adaptive_order(alg::ImplicitTaylor{P, Q}) where {P, Q}
    if is_mu_taylor(alg)
        return 0 # haven't implemented embedded method for μ-Taylor yet
    else
        return P + Q - 1
    end
end

function normalized_pade(p::Int, q::Int)
    RI = Rational{Int}
    N = p + q
    c = Vector{RI}(undef, N + 1)
    c[1] = 1 // 1
    for k in 1:N
        c[k + 1] = c[k] // k
    end
    if q == 0
        dres = RI[]
    else
        A = zeros(RI, q, q)
        rhs = zeros(RI, q)
        for i in 1:q
            k = p + i
            rhs[i] = -c[k + 1]
            for j in 1:q
                idx = k - j
                A[i, j] = idx >= 0 ? c[idx + 1] : 0 // 1
            end
        end
        dres = A \ rhs
    end
    n = zeros(RI, p + 1)
    for k in 0:p
        n[k + 1] = c[k + 1]
        for j in 1:min(k, q)
            n[k + 1] += dres[j] * c[k - j + 1]
        end
    end
    d = vcat(1 // 1, dres)
    normalized_n = [x * factorial(k - 1) for (k, x) in enumerate(n)]
    normalized_d = [x * factorial(k - 1) for (k, x) in enumerate(d)]
    return normalized_n, normalized_d
end
