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

# Error estimate (see implicit_taylor_error!):
#   * order ≥ 3 Taylor-Padé: filtered residual W⁻¹ r1 of an embedded formula of
#     one order less, W = (I - hJ/(P+Q))^Q. For Q = 2 the embedded formula is
#     (P - 1, Q) ((0, 1) for (0, 2)): both implicit polynomials have the same
#     degree, so the filtered Newton error is O(δ). For Q ≥ 3 it is (P, Q - 1):
#     a formula with the same implicit degree shares the spurious roots that
#     approach the zeros of c_Q for large h, and its residual vanishes there too.
#   * μ-Taylor: the degree-P term at the expansion center.
#   * second order ((1, 1), (0, 2)): third divided differences of the last four
#     states; first order: second divided differences of the last three.
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

# The Newton factorization is kept while no row of the newly assembled matrix
# differs from the factorized one by more than this, relative to
# max(1, row sum of magnitudes) (see newton_matrix_outdated).
matrix_cutoff(::ImplicitTaylor) = 1 // 20

# Same for the error filter; the estimate only needs to be right within a factor.
filter_cutoff(::ImplicitTaylor) = 1 // 5

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

# Weights of the degree-P Taylor increment: Σ_{k=1}^{P} c_k h^k.
taylor_increment_coefficients(::Val{P}) where {P} = ntuple(k -> k == 1 ? 0.0 : 1.0, P + 1)

# Weights selecting the single highest term c_P h^P of a degree-P polynomial.
taylor_last_coefficients(::Val{P}) where {P} = ntuple(k -> k == P + 1 ? 1.0 : 0.0, P + 1)
