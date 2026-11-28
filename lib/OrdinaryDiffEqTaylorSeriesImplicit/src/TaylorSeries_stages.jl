# In derivative-stage form the scaled coefficients v_k = h^k c_k(u) are unknowns:
#   R_0 = z + Σ_k b_k v_k - A_inc(uprev),                         u = uprev + z,
#   R_k = v_k - (h/k) F_{k-1}(u, v_1, …, v_{k-1}),               k = 1:Q.
# F_j is the j-th coefficient of f along the stage jet. Simplified Newton
# freezes J = f_u(u,t+h) and eliminates z, giving (I + C ⊗ hJ) δv = r,
# where r_1 = R_1 + hJ R_0, r_k = R_k for k ≥ 2, and δz = R_0 - Σ_k b_k δv_k.
# C has first row b_1 … b_Q and subdiagonal -1/k. Its eigendecomposition splits
# the system into I + λ_i hJ factors, one real per real eigenvalue and one
# complex per conjugate pair.

# C = T Λ T⁻¹ for the weights b_k (normalized Padé denominator coefficients), split into real
# eigenvalues and one representative (imag > 0) of each conjugate pair. The eigenvalues are
# -1/r over the roots r of D; the subdiagonal gives the eigenvectors in closed form,
# v_1 = 1, v_k = -v_{k-1}/(kλ).
const StageEigensystem = Tuple{
    Vector{Float64}, Vector{ComplexF64}, Matrix{Float64},
    Matrix{ComplexF64}, Matrix{Float64}, Matrix{ComplexF64},
}
# The eigensystem depends on the method alone: computed once per set of weights.
const STAGE_EIGENSYSTEMS = Dict{Vector{Float64}, StageEigensystem}()
const STAGE_EIGENSYSTEMS_LOCK = ReentrantLock()
function stage_eigensystem(weights)
    key = Float64[w for w in weights]
    return lock(STAGE_EIGENSYSTEMS_LOCK) do
        get!(() -> compute_stage_eigensystem(weights), STAGE_EIGENSYSTEMS, key)
    end::StageEigensystem
end
function compute_stage_eigensystem(weights)
    Q = length(weights) - 1
    roots = polynomial_roots([Float64(weights[k + 1]) / factorial(k) for k in 1:Q])
    λ = [-1 / r for r in roots]
    T = zeros(ComplexF64, Q, Q)
    for (i, l) in enumerate(λ)
        T[1, i] = 1
        for k in 2:Q
            T[k, i] = -T[k - 1, i] / (k * l)
        end
    end
    Tinv = small_inverse(T)
    isreal_λ = [abs(imag(l)) <= 1.0e-10 * abs(l) for l in λ]
    real_idx = findall(isreal_λ)
    cplx_idx = findall(i -> !isreal_λ[i] && imag(λ[i]) > 0, eachindex(λ))
    length(real_idx) + 2length(cplx_idx) == Q || error("stage eigensystem: unexpected eigenvalues $λ")
    return (
        real.(λ[real_idx]), λ[cplx_idx],
        real.(T[:, real_idx]), T[:, cplx_idx], real.(Tinv[real_idx, :]), Tinv[cplx_idx, :],
    )
end

# Inverse of a small matrix by Gauss-Jordan elimination with partial pivoting.
function small_inverse(A)
    Q = size(A, 1)
    M = hcat(copy(A), [i == j ? one(eltype(A)) : zero(eltype(A)) for i in 1:Q, j in 1:Q])
    for c in 1:Q
        p = argmax(i -> abs(M[i, c]), c:Q)
        if p != c
            for j in axes(M, 2)
                M[c, j], M[p, j] = M[p, j], M[c, j]
            end
        end
        pivot = M[c, c]
        for j in axes(M, 2)
            M[c, j] /= pivot
        end
        for i in 1:Q
            i == c && continue
            factor = M[i, c]
            for j in axes(M, 2)
                M[i, j] -= factor * M[c, j]
            end
        end
    end
    return M[:, (Q + 1):end]
end

# The factors (I + λ_i hJ) of the simplified Newton matrix, their solvers and work vectors.
struct StageFactors{T, LR, LC}
    λr::Vector{Float64}
    λc::Vector{ComplexF64}
    Tr::Matrix{Float64}        # Q × (real roots): columns of T
    Tc::Matrix{ComplexF64}     # Q × (pairs)
    Tir::Matrix{Float64}       # rows of T⁻¹
    Tic::Matrix{ComplexF64}
    Ar::Vector{Matrix{T}}
    Ac::Vector{Matrix{Complex{T}}}
    lr::Vector{LR}
    lc::Vector{LC}
    wr::Vector{Vector{T}}      # right-hand sides / solutions
    wc::Vector{Vector{Complex{T}}}
    direct::Bool               # factorize and solve in place with dense_lu!/dense_solve!
    pivr::Vector{Vector{Int}}  # their pivots
    pivc::Vector{Vector{Int}}
    h::Base.RefValue{T}        # distance of the current factorizations
    hJ::Matrix{T}              # h J of the current factorizations
    rowchange::Vector{T}       # row sums of |hJ_new - hJ| and of the off-diagonal |hJ_new|
    rowsum::Vector{T}
end

# Use the in-place LU kernels for small dense factors with the default solver.
const DIRECT_LU_MAX = 24

function build_stage_factors(
        weights, J::Matrix{T}, linsolve_alg, linear_u, p, t;
        direct = linsolve_alg === nothing && size(J, 1) <= DIRECT_LU_MAX
    ) where {T}
    λr, λc, Tr, Tc, Tir, Tic = stage_eigensystem(weights)
    n = size(J, 1)
    make(A, b) = init(
        LinearProblem(A, b, (nothing, linear_u, p, t); u0 = copy(b)), linsolve_alg,
        alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        assumptions = LinearSolve.OperatorAssumptions(true)
    )
    Ar = [zeros(T, n, n) for _ in λr]
    Ac = [zeros(Complex{T}, n, n) for _ in λc]
    wr = [zeros(T, n) for _ in λr]
    wc = [zeros(Complex{T}, n) for _ in λc]
    lr = [make(Ar[i], wr[i]) for i in eachindex(λr)]
    lc = [make(Ac[i], wc[i]) for i in eachindex(λc)]
    pivr = [zeros(Int, n) for _ in λr]
    pivc = [zeros(Int, n) for _ in λc]
    return StageFactors(
        λr, λc, Tr, Tc, Tir, Tic, Ar, Ac, lr, lc, wr, wc,
        direct, pivr, pivc,
        Ref(convert(T, NaN)), zeros(T, n, n), zeros(T, n), zeros(T, n)
    )
end

# Empty factors for a cache without a derivative-stage system.
function empty_stage_factors(::Type{T}, linsolve_alg, linear_u, p, t) where {T}
    make(A, b) = init(
        LinearProblem(A, b, (nothing, linear_u, p, t); u0 = copy(b)), linsolve_alg,
        alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        assumptions = LinearSolve.OperatorAssumptions(true)
    )
    Ar, wr = zeros(T, 1, 1), zeros(T, 1)
    Ac, wc = zeros(Complex{T}, 1, 1), zeros(Complex{T}, 1)
    lr, lc = [make(Ar, wr)], [make(Ac, wc)]
    return StageFactors(
        Float64[], ComplexF64[], zeros(0, 0), zeros(ComplexF64, 0, 0),
        zeros(0, 0), zeros(ComplexF64, 0, 0), [Ar], [Ac], lr, lc, [wr], [wc],
        false, [Int[]], [Int[]], Ref(convert(T, NaN)), zeros(T, 0, 0), T[], T[]
    )
end

# A partial LU has incomplete pivots and must never reach back substitution.
function dense_lu!(A::AbstractMatrix, piv::AbstractVector{Int})
    lu_factor!(A, piv) ||
        throw(TaylorLinearSolveFailure("ImplicitTaylor encountered a singular or nonfinite LU factor"))
    return A
end
dense_solve!(A::AbstractMatrix, piv::AbstractVector{Int}, b::AbstractVector) =
    check_finite_solution(lu_solve!(A, piv, b))

function dense_solve!(A::AbstractMatrix, piv::AbstractVector{Int}, B::AbstractMatrix)
    for j in axes(B, 2)
        dense_solve!(A, piv, @view(B[:, j]))
    end
    return B
end

nfactorizations(S::StageFactors) = length(S.λr) + length(S.λc)

# (I + λ_i hJ) for every root; the factorization happens on the next solve.
function form_stage_factors!(S::StageFactors, J, h)
    n = size(J, 1)
    for (i, λ) in enumerate(S.λr)
        A = S.Ar[i]
        s = λ * h
        @inbounds for j in 1:n, k in 1:n
            A[k, j] = s * J[k, j]
        end
        @inbounds for k in 1:n
            A[k, k] += one(eltype(A))
        end
    end
    for (i, λ) in enumerate(S.λc)
        A = S.Ac[i]
        s = λ * h
        @inbounds for j in 1:n, k in 1:n
            A[k, j] = s * J[k, j]
        end
        @inbounds for k in 1:n
            A[k, k] += one(eltype(A))
        end
    end
    S.h[] = h
    S.hJ .= h .* J
    return nothing
end

# Check the row-relative change of every factor I + λ_i hJ. The off-diagonal
# differences are |λ_i| times the corresponding differences in hJ.
function stage_factors_outdated(S::StageFactors, J, h, cutoff)
    isnan(S.h[]) && return true
    n = size(J, 1)
    change, offdiag = S.rowchange, S.rowsum
    fill!(change, zero(eltype(change)))
    fill!(offdiag, zero(eltype(offdiag)))
    @inbounds for j in 1:n
        for k in 1:n
            x = h * J[k, j]
            change[k] += abs(x - S.hJ[k, j])
            k == j || (offdiag[k] += abs(x))
        end
    end
    @inbounds for k in 1:n
        x = h * J[k, k]
        for λ in S.λr
            abs(λ) * change[k] > cutoff * max(one(x), abs(1 + λ * x) + abs(λ) * offdiag[k]) && return true
        end
        for λ in S.λc
            abs(λ) * change[k] > cutoff * max(one(x), abs(1 + λ * x) + abs(λ) * offdiag[k]) && return true
        end
    end
    return false
end

@inline function solve_stage_factor!(integrator, linsolve, A, w, factorize)
    # The linear solver's solution buffer must not alias its right-hand side.
    result = factorize ? checked_dolinsolve(integrator, linsolve; A, b = w) :
        checked_dolinsolve(integrator, linsolve; b = w)
    return copyto!(w, result.u)
end

# dv ← (I + C ⊗ hJ)⁻¹ rb, factorizing the formed factors first when `factorize`.
function solve_stage_factors!(integrator, S::StageFactors, dv, rb, n, factorize)
    Q = size(S.Tr, 1)
    for i in eachindex(S.λr)
        w = S.wr[i]
        @inbounds for k in 1:n
            acc = zero(eltype(w))
            for m in 1:Q
                acc += S.Tir[i, m] * rb[(m - 1) * n + k]
            end
            w[k] = acc
        end
        if S.direct
            factorize && dense_lu!(S.Ar[i], S.pivr[i])
            dense_solve!(S.Ar[i], S.pivr[i], w)
        else
            solve_stage_factor!(integrator, S.lr[i], S.Ar[i], w, factorize)
        end
    end
    for i in eachindex(S.λc)
        w = S.wc[i]
        @inbounds for k in 1:n
            acc = zero(eltype(w))
            for m in 1:Q
                acc += S.Tic[i, m] * rb[(m - 1) * n + k]
            end
            w[k] = acc
        end
        if S.direct
            factorize && dense_lu!(S.Ac[i], S.pivc[i])
            dense_solve!(S.Ac[i], S.pivc[i], w)
        else
            solve_stage_factor!(integrator, S.lc[i], S.Ac[i], w, factorize)
        end
    end
    @inbounds for m in 1:Q, k in 1:n
        acc = zero(eltype(dv))
        for i in eachindex(S.λr)
            acc += S.Tr[m, i] * S.wr[i][k]
        end
        for i in eachindex(S.λc)
            acc += 2 * real(S.Tc[m, i] * S.wc[i][k])
        end
        dv[(m - 1) * n + k] = acc
    end
    return dv
end

# R_0 into r0 and the right-hand side r_1 = R_1 + hJ R_0, r_k = R_k (k ≥ 2) into rb. The
# elimination term uses the matrix of the factors, `scale * J` (h0 J0 for kept factors).
function stage_residual!(rb, r0, z, V, Fst, J, A, weights, h, scale = h)
    n, Q = size(V, 1), size(V, 2)
    @inbounds for i in 1:n
        acc = z[i] - A[i]
        for k in 1:Q
            acc += weights[k + 1] * V[i, k]
        end
        r0[i] = acc
    end
    @inbounds for k in 1:Q
        s = h / k
        b0 = (k - 1) * n
        for i in 1:n
            acc = V[i, k] - s * Fst[i, k]
            if k == 1
                for j in 1:n
                    acc += scale * J[i, j] * r0[j]
                end
            end
            rb[b0 + i] = acc
        end
    end
    return rb
end

# Stage predictor: the previous step's converged stages v_k = h0^k c_k(uprev) rescaled to h.
function stage_predictor!(V, V0, valid, h0, h)
    if valid && !iszero(h0)
        ratio = h / h0
        scale = one(ratio)
        @inbounds for k in axes(V, 2)
            scale *= ratio
            for i in axes(V, 1)
                V[i, k] = scale * V0[i, k]
            end
        end
    else
        fill!(V, zero(eltype(V)))
    end
    return V
end

# δz = R_0 - Σ_k b_k δv_k from the stage corrections.
function stage_increment_correction!(dz, r0, dv, V, weights)
    n, Q = size(V, 1), size(V, 2)
    @inbounds for i in 1:n
        acc = r0[i]
        for k in 1:Q
            acc -= weights[k + 1] * dv[(k - 1) * n + i]
        end
        dz[i] = acc
    end
    return dz
end

# Scale the stage-value correction by the integration tolerances.
stage_correction_norm!(atmp, dz, uprev, u, abstol, reltol, internalnorm, t) =
    (
    @.. broadcast = false atmp = calculate_residuals(abs(dz), uprev, u, abstol, reltol, internalnorm, t);
    internalnorm(atmp, t)
)

# Recompute the equations at the updated iterate. A small correction to z can
# hide large, cancelling corrections to the auxiliary Taylor coefficients.
function check_derivative_defect!(integrator, cache, t, tc, h, uprev, weights)
    (; κ, uintermediate, utilde, V, Fst, tmp, tmp_embedded, stage_error, Mtmp) = cache
    (; adaptive, abstol, reltol, internalnorm) = integrator.opts
    Q = size(V, 2)
    all(isfinite, V) || return oftype(κ, Inf), false
    cache.stage_values(Fst, uintermediate, V, tc, h)
    all(isfinite, Fst) || return oftype(κ, Inf), false

    filter_ready = false
    paired = false
    try_reduced = adaptive && (iszero(cache.reduced_check_h) || 2 * abs(h) < cache.reduced_check_h)
    if adaptive
        cache.polynomial_check(tmp, tmp_embedded, cache.Jcheck, uintermediate, tc, h)
        integrator.stats.njacs += 1
        cache.ncheckjac += 1
        # Reuse the filter for both the reduced equation and the embedded
        # residual when the backend supports two right-hand sides.
        if try_reduced
            @. stage_error = utilde + tmp - cache.rhs
            if all(isfinite, stage_error)
                paired = cache.linsolve_filter.b isa AbstractMatrix
                if paired
                    @. tmp_embedded = utilde + tmp_embedded - cache.rhs_embedded
                end
                filter_error!(
                    integrator, cache, stage_error, Q, weights[2], h, Mtmp,
                    paired ? tmp_embedded : nothing
                )
                filter_ready = true
                reduced = stage_correction_norm!(
                    cache.atmp, stage_error, uprev, uintermediate, abstol, reltol, internalnorm, t
                )
                if isfinite(reduced) && reduced <= κ
                    cache.reduced_check_h = zero(h)
                    return reduced, paired
                end
            end
        end
    end

    # Validate each lifted equation separately. R_0 is already in solution
    # units. R_k contributes to the eliminated residual through degree Q-k,
    # so apply k filter factors for adaptive steps. A fixed step needs no LTE
    # filter, which can have poles even when the Padé equation is invertible.
    r = cache.r0
    @inbounds for i in eachindex(r)
        value = utilde[i] - cache.rhs[i]
        for k in 1:Q
            value += weights[k + 1] * V[i, k]
        end
        r[i] = value
    end
    all(isfinite, r) || return oftype(κ, Inf), false
    defect = stage_correction_norm!(cache.atmp, r, uprev, uintermediate, abstol, reltol, internalnorm, t)
    for k in 1:Q
        @inbounds for i in eachindex(r)
            r[i] = V[i, k] - (h / k) * Fst[i, k]
        end
        all(isfinite, r) || return oftype(κ, Inf), false
        if adaptive
            if filter_ready
                apply_filter_factors!(integrator, cache, r, k)
            else
                filter_error!(integrator, cache, r, Q, weights[2], h, Mtmp; count = k)
                filter_ready = true
            end
        end
        defect = max(defect, stage_correction_norm!(cache.atmp, r, uprev, uintermediate, abstol, reltol, internalnorm, t))
    end
    # The reduced jet can be rounding-limited at large distances. Reconsider
    # it after halving the distance; the lifted check remains mandatory here.
    if try_reduced && isfinite(defect) && defect <= κ
        cache.reduced_check_h = abs(h)
    end
    return defect, paired
end

function apply_stage_correction!(z, V, dz, dv)
    n, Q = size(V, 1), size(V, 2)
    @inbounds for i in 1:n
        z[i] -= dz[i]
    end
    @inbounds for k in 1:Q, i in 1:n
        V[i, k] -= dv[(k - 1) * n + i]
    end
    return nothing
end

# :auto retries a failed reduced solve at the same step size in derivative form.
# In adaptive mode it also switches when slow convergence prevents proposed
# growth, and probes the reduced form again after fast steps, with backoff.
function hand_over_to_stages!(cache, adaptive)
    cache.recovery_since_accept = cache.stage_recovery
    if adaptive
        cache.stage_trial && (cache.stage_probe *= 2)
        cache.use_stages = true
        cache.stage_handover = true
    end
    cache.stage_trial = false
    cache.stage_streak = 0
    return nothing
end

# Mode bookkeeping at an accepted step (adaptive :auto), with the controller's
# effective q before Newton feedback (q < 1: growth).
function stage_mode_feedback!(state, q)
    held = state.newton_limit || state.recovery_since_accept
    if state.use_stages
        if state.newton_failed || held
            state.stage_streak = 0
        else
            state.stage_streak += 1
            if state.stage_streak >= state.stage_probe
                state.use_stages = false
                state.stage_trial = true
                state.stage_streak = 0
            end
        end
    elseif held && q < one(q)
        state.stage_trial && (state.stage_probe *= 2)
        state.use_stages = true
        state.stage_trial = false
        state.stage_streak = 0
    else
        state.stage_trial && (state.stage_probe = 1)
        state.stage_trial = false
    end
    state.stage_handover = false
    return nothing
end

# Estimate how a stage correction δ would leak into the next explicit Padé
# polynomial. Approximate the next Jacobian by J + (J - Jprev) and use
# Σ_k a_k h^k/k! [(2J - Jprev)^k - J^k] δ, then apply one error-filter factor
# and scale by the integration tolerances. The leak vanishes for constant J.
# The stopping rule scales the measured leak by η to estimate the remaining error.
function stage_leak_norm!(
        integrator, cache, δ, J, Jprev, h, Q, explicit_weights, b1, uprev, u,
        abstol, reltol, internalnorm, t
    )
    (; leak_p, leak_q, leak_w, stage_error) = cache
    P = length(explicit_weights) - 1
    copyto!(leak_p, δ)              # (J + ΔJ)^k δ,  J + ΔJ = 2J - Jprev
    copyto!(leak_q, δ)              # J^k δ
    fill!(stage_error, zero(eltype(stage_error)))
    hk = one(h)
    kfact = 1
    for k in 1:P
        hk *= h
        kfact *= k
        c = explicit_weights[k + 1] * hk / kfact
        matvec!(leak_w, J, leak_q, Jprev, zero(c))
        copyto!(leak_q, leak_w)
        matvec!(leak_w, J, leak_p, Jprev, one(c))    # (2J - Jprev) p = J p + (J p - Jprev p)
        copyto!(leak_p, leak_w)
        @. stage_error += c * (leak_p - leak_q)
    end
    # Project the fast part out with one factor of the error filter (factorized here and
    # reused by the error estimate).
    filter_error!(integrator, cache, stage_error, Q, b1, h, J; count = 1)
    return stage_correction_norm!(cache.atmp, stage_error, uprev, u, abstol, reltol, internalnorm, t)
end

# w ← J x + α (J - B) x, column by column (α = 0: J x; α = 1: (2J - B) x).
function matvec!(w, J, x, B, α)
    n = length(x)
    fill!(w, zero(eltype(w)))
    @inbounds for j in 1:n
        xj = x[j]
        if iszero(α)
            @simd for i in 1:n
                w[i] += J[i, j] * xj
            end
        else
            @simd for i in 1:n
                w[i] += (J[i, j] + α * (J[i, j] - B[i, j])) * xj
            end
        end
    end
    return w
end

# Solve the derivative-stage equations with f_u evaluated at the predictor.
# Adaptive solves reuse factors within the row cutoff and refresh when needed.
function derivative_stage!(integrator, cache, t, dt, uprev, u, tc, distance, weights, new_step)
    (; κ, uintermediate, rhs, utilde) = cache
    (; V, V0, Fst, Jst, rb, r0, dz, stage_factors) = cache
    alg = unwrap_alg(integrator, true)
    (; internalnorm, abstol, reltol, adaptive) = integrator.opts
    (; maxiters) = alg
    stats = integrator.stats
    h = distance
    θ_slow = convert(typeof(κ), 1 // 2)
    κ_second = κ * convert(typeof(κ), 1 // 10)
    tiny = κ * convert(typeof(κ), 1 // 1000)
    ηleak = leak_stop_rate(cache, h)   # leak rate estimate: an earlier step's until measured here
    leak_measured = false
    dv = cache.dv
    n = length(uprev)
    J = @view Jst[:, :, 1]
    cache.embedded_filtered = false
    # A retry at the same predictor can reuse its Jacobian.
    evaluate = !(cache.stage_jac_valid && (cache.constant_stage_jacobian || (adaptive && !new_step)))
    force = !adaptive && !cache.constant_stage_jacobian
    η = ndw = ndwprev = θ = one(κ)
    converged = false
    measured = false
    cache_retry = false
    iter = 0
    total_iter = 0
    # Cache the start-of-step Jacobian for the leak estimate.
    P, Q = _unwrap_val(alg.order), _unwrap_val(alg.order_q)
    leak_test = P >= 1 && !cache.constant_stage_jacobian
    if leak_test && (new_step || !cache.jprev_valid)
        copyto!(uintermediate, uprev)
        cache.stage_jac(Fst, cache.Jprev, uintermediate, V, t, h)
        stats.njacs += 1
        cache.jprev_valid = true
    end
    explicit_weights = pade_increment_coefficients(Val(P), Val(Q))[1]
    Jprev = @view cache.Jprev[:, :, 1]
    for attempt in 1:2
        taylor_initial_guess!(integrator, cache, alg, t, dt, uprev, true)
        @. uintermediate = uprev + utilde
        stage_predictor!(V, V0, cache.stages_valid, cache.stages_h, h)
        cache.status = VerySlowConvergence
        η = one(η)
        iter = 0
        jac_iter = 1
        θ_from = 3
        refresh = false       # refresh the matrix at the current iterate before the next solve
        refreshed = false     # at most one refresh per attempt
        fresh_matrix = false  # the matrix in use was evaluated in this attempt
        leak_prev = zero(κ)   # leak of the last tested correction and its iteration
        leak_iter = 0
        validation_previous = oftype(κ, Inf)
        while iter < maxiters
            iter += 1
            total_iter += 1
            stats.nnonliniter += 1
            factorize = false
            if iter == 1 || refresh || !adaptive
                jac_now = evaluate || refresh || (!adaptive && !cache.constant_stage_jacobian)
                if jac_now
                    cache.stage_jac(Fst, Jst, uintermediate, V, tc, h)
                    stats.njacs += 1
                    cache.ncoefjac += 1
                    cache.stage_jac_valid = true
                else
                    cache.stage_values(Fst, uintermediate, V, tc, h)
                end
                outdated = !cache.lu_stages_valid ||
                    (
                    adaptive ? stage_factors_outdated(stage_factors, J, h, matrix_cutoff(alg)) :
                        stage_factors.h[] != h
                )
                if force || refresh || outdated
                    form_stage_factors!(stage_factors, J, h)
                    cache.lu_stages_valid = true
                    factorize = true
                    # Kept factors are not fresh: if they fail, the attempt is retried.
                    jac_now && (fresh_matrix = true)
                end
                if refresh
                    cache.nendpoint_refresh += 1
                    cache.recovery_since_accept = true
                    refresh = false
                end
                jac_iter = iter
                θ_from = iter + 2
            else
                cache.stage_values(Fst, uintermediate, V, tc, h)
            end
            stage_residual!(rb, r0, utilde, V, Fst, stage_factors.hJ, rhs, weights, h, one(h))
            factorize && (stats.nw += nfactorizations(stage_factors))
            solve_stage_factors!(integrator, stage_factors, dv, rb, n, factorize)
            stats.nsolve += nfactorizations(stage_factors)
            stage_increment_correction!(dz, r0, dv, V, weights)
            iter > 1 && (ndwprev = ndw)
            ndw = stage_correction_norm!(cache.atmp, dz, uprev, uintermediate, abstol, reltol, internalnorm, t)
            if !isfinite(ndw)
                cache.status = Divergence
                break
            end
            apply_stage_correction!(utilde, V, dz, dv)
            @. uintermediate = uprev + utilde

            bad_convergence = false
            slow = false
            if iter >= θ_from
                θ = ndw / ndwprev
                bad_convergence = θ > 2 || (adaptive && θ >= one(θ))
                slow = adaptive && θ > θ_slow
                η = θ / (1 - θ)
            end
            candidate = ndw <= tiny || (adaptive ? iter == jac_iter + 1 : iter >= 2) && ndw < κ_second ||
                (iter >= θ_from && η >= zero(η) && η * ndw < κ)
            if candidate && !bad_convergence && leak_test
                # Bound the remaining correction’s effect on the next explicit polynomial.
                leak = stage_leak_norm!(
                    integrator, cache, dz, J, Jprev, h, Q, explicit_weights,
                    weights[2], uprev, uintermediate, abstol, reltol, internalnorm, t
                )
                if leak_iter > 0 && leak_prev > zero(leak_prev) && isfinite(leak)
                    # Per-iteration contraction of the leak since the last tested correction.
                    θf = (leak / leak_prev)^(1 // (iter - leak_iter))
                    ηleak = θf < θ_slow ? θf / (1 - θf) : one(θf)
                    leak_measured = true
                end
                leak_prev, leak_iter = leak, iter
                candidate = isfinite(leak) && ηleak * leak < κ
                candidate || (cache.nleak += 1)
            end
            if candidate && !bad_convergence
                defect, paired = check_derivative_defect!(integrator, cache, t, tc, h, uprev, weights)
                if isfinite(defect) && defect <= κ
                    cache.embedded_filtered = paired
                    converged = true
                    break
                end
                cache.nvalidation_fail += 1
                if !isfinite(defect)
                    cache.status = Divergence
                    break
                end
                # Let the verified residual contract and refresh once if it stalls.
                η = one(η)
                θ_from = iter + 2
                if !refreshed && defect >= θ_slow * validation_previous
                    refreshed = refresh = true
                end
                validation_previous = defect
                continue
            end
            if bad_convergence
                cache.status = Divergence
                break
            end
            if slow
                # A kept or stale matrix: refresh it at the current iterate once.
                if !refreshed
                    refreshed = true
                    refresh = true
                    continue
                end
                cache.status = VerySlowConvergence
                break
            end
            # The remaining budget cannot reach the tolerance at this θ.
            if adaptive && iter >= θ_from && θ < one(θ) &&
                    ndw * θ^(maxiters - iter) > κ * (1 - θ)
                cache.status = VerySlowConvergence
                break
            end
        end
        measured = iter >= θ_from
        # Retry once from the predictor with a fresh matrix unless this attempt already had one.
        (converged || fresh_matrix || attempt == 2) && break
        cache.nrestart += 1
        cache_retry = true
        evaluate = true
        force = true
    end
    cache.iter = total_iter
    if !converged
        cache.stage_jac_valid = false
        cache.lu_stages_valid = false
        cache.newton_failed = true
        adaptive && record_newton_failure!(cache, dt)
        integrator.force_stepfail = true
        stats.nnonlinconvfail += 1
        return false
    end
    cache.status = Convergence
    fast = !measured || zero(η) <= η <= fast_convergence_cutoff(alg)
    cache_retry && !fast && (cache.recovery_since_accept = true)
    leak_measured && (cache.ηleak = ηleak; cache.ηleak_h = abs(h))
    record_newton_convergence!(cache, measured, η, alg)
    # Only validated candidates contribute predictor history. Validation also
    # prepared the embedded polynomial and the endpoint Jacobian for adaptivity.
    copyto!(V0, V)
    cache.stages_h = h
    cache.stages_valid = true
    if !adaptive
        copyto!(cache.Mtmp, @view(Jst[:, :, 1]))
    end
    u .= uintermediate
    return true
end
