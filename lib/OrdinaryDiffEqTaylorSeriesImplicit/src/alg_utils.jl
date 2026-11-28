alg_extrapolates(::ImplicitTaylor) = true

# Required by the algorithm interface; the Taylor kernels do their own differentiation.
alg_autodiff(::ImplicitTaylor) = AutoForwardDiff()

is_mu_taylor(::ImplicitTaylor{P, Q}) where {P, Q} = Q == 0

is_taylor_pade(::ImplicitTaylor{P, Q}) where {P, Q} = Q != 0

function alg_order(alg::ImplicitTaylor{P, Q}) where {P, Q}
    if is_mu_taylor(alg) # μ-Taylor method
        μ = alg.μ
        # Leading error ∝ a - b; at the next power of h two elementary
        # differentials appear, with coefficients a and c.
        a = (-μ)^(P + 1)
        b = (1 - μ)^(P + 1)
        leading = a - b
        c = ((-μ)^(P + 2) - (1 - μ)^(P + 2)) / (P + 2) + μ * leading
        project(x) = alg.real_function ? real(x) : x
        tol = 8 * eps(real(float(one(μ))))
        cancelled(x, scale) = abs(project(x)) <= tol * scale
        cancelled(leading, abs(a) + abs(b)) || return P
        next_scale = abs(a) + abs(b) + abs(c)
        return cancelled(a, next_scale) && cancelled(c, next_scale) ? P + 2 : P + 1
    else # Taylor-Padé method
        return P + Q
    end
end

# Higher-order Padé methods use an embedded formula one order lower. For Q ≥ 3
# reduce the implicit degree so the estimator can detect spurious roots of the
# main formula. Low orders use divided differences and μ-Taylor uses its last term.
uses_embedded(::ImplicitTaylor{P, Q}) where {P, Q} = Q != 0 && P + Q >= 3
uses_divided3(::ImplicitTaylor{P, Q}) where {P, Q} = Q != 0 && P + Q == 2
uses_last_term(::ImplicitTaylor{P, Q}) where {P, Q} = Q == 0 && P >= 2

function alg_adaptive_order(alg::ImplicitTaylor{P, Q}) where {P, Q}
    if uses_embedded(alg)
        return P + Q - 1
    elseif uses_divided3(alg)
        return 2
    elseif uses_last_term(alg)
        return P - 1
    else
        return 1
    end
end

# Padé error constant: local error C h^{P+Q+1} y^{(P+Q+1)} on y' = λy.
pade_error_constant(P, Q) = factorial(P) * factorial(Q) // (factorial(P + Q) * factorial(P + Q + 1))

# Keep dt constant when the controller would change it only slightly.
OrdinaryDiffEqCore.qsteady_max_default(::ImplicitTaylor) = 6 // 5

function OrdinaryDiffEqCore.default_controller(QT, alg::ImplicitTaylor)
    return alg.controller === :Predictive ? OrdinaryDiffEqCore.PredictiveController(QT, alg) :
        OrdinaryDiffEqCore.PIController(QT, alg)
end
# Plain safety factor γ: the Newton feedback replaces the iteration-count factor.
OrdinaryDiffEqCore.fac_default_gamma(::ImplicitTaylor) = true

# η = θ/(1 - θ) below this counts as fast convergence; dt grows only after it.
fast_convergence_cutoff(::ImplicitTaylor) = 1 // 5

# Relative row-change cutoff for reusing a Newton factorization.
matrix_cutoff(::ImplicitTaylor) = 1 // 20

# Relative row-change cutoff for reusing the error filter.
filter_cutoff(::ImplicitTaylor) = 1 // 5

# If the embedded implicit degree Q₁ exceeds P, reusing the filter can amplify
# stiff Newton error by s^(Q₁-P), where s is the scaled Jacobian row norm.
# Return (A, m) so the absolute row change is capped at A s^(1-m) for s > 1.
function filter_leak_limit(alg::ImplicitTaylor{P, Q}, κ) where {P, Q}
    growth, m = filter_leak_growth(Val(P), Val(Q))
    m > 0 || return (zero(κ), 0)
    return (convert(typeof(κ), filter_cutoff(alg)) / (κ * convert(typeof(κ), growth)), m)
end


embedded_degrees(P, Q) = Q >= 3 ? (P, Q - 1) : (P >= 1 ? (P - 1, Q) : (0, Q - 1))

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

# Padé coefficients, computed once per (P, Q) at compile time; constant terms are
# dropped because both polynomials are used in increment form.
@generated function pade_increment_coefficients(::Val{P}, ::Val{Q}) where {P, Q}
    a, b = normalized_pade(P, Q)
    a[1] = 0
    b[1] = 0
    return QuoteNode((Tuple(Float64.(a)), Tuple(Float64.(b))))
end

# (b'_{Q₁}/a_P) (Q/|b₁|)^m and m = Q₁ - P, at compile time; m = 0 when no embedded
# formula is used or its degree does not exceed P.
@generated function filter_leak_growth(::Val{P}, ::Val{Q}) where {P, Q}
    (Q != 0 && P + Q >= 3) || return QuoteNode((0.0, 0))
    P1, Q1 = embedded_degrees(P, Q)
    m = Q1 - P
    m > 0 || return QuoteNode((0.0, 0))
    a, b = normalized_pade(P, Q)
    b1 = normalized_pade(P1, Q1)[2]
    aP = P == 0 ? 1.0 : Float64(a[P + 1])
    growth = (Q / abs(Float64(b[2])))^m * abs(Float64(b1[Q1 + 1])) / aP
    return QuoteNode((growth, m))
end

# Weights of the degree-P Taylor increment: Σ_{k=1}^{P} c_k h^k.
taylor_increment_coefficients(::Val{P}) where {P} = ntuple(k -> k == 1 ? 0.0 : 1.0, P + 1)

# Weights selecting the single highest term c_P h^P of a degree-P polynomial.
taylor_last_coefficients(::Val{P}) where {P} = ntuple(k -> k == P + 1 ? 1.0 : 0.0, P + 1)
