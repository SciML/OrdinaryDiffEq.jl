# Use D(hJ) = Π_i (I - hJ/r_i) as a simplified Newton matrix, with
# J = ∂f/∂u and r_i the roots of the Padé denominator. Each factor avoids
# forming (hJ)^Q. Higher derivatives of f still enter the nonlinear residual.
# Available for adaptive Taylor-Padé methods with Q ≥ 3 on real dense problems.

struct FactoredNewton{T}
    real_roots::Vector{T}
    complex_roots::Vector{Complex{T}}      # one root of each conjugate pair (imag > 0)
    real_factors::Vector{Matrix{T}}        # LU factors in place
    real_pivots::Vector{Vector{Int}}
    complex_factors::Vector{Matrix{Complex{T}}}
    complex_pivots::Vector{Vector{Int}}
    cbuf::Vector{Complex{T}}
    correction::Vector{T}
    distance::Base.RefValue{T}             # distance of the current factors (NaN: none)
    enabled::Bool                          # newton_matrix = :auto or :factored (else inert, size 0)
    auto::Bool                             # newton_matrix = :auto
    active::Base.RefValue{Bool}            # the factored iteration is in use
    rounding::Base.RefValue{T}             # eps·Π cond(factor) of the current factors
    rate::Base.RefValue{T}                 # last measured contraction rate of the factored iteration
    # Work vectors of the condition estimates.
    rx::Vector{T}
    ry::Vector{T}
    cx::Vector{Complex{T}}
    cy::Vector{Complex{T}}
    # The assembled iteration's dt cap and cap probing while the factored one runs.
    assembled_cap::Base.RefValue{T}
    assembled_probe::Base.RefValue{NTuple{3, Int}}
end

# Roots of 1 + Σ_{k=1}^Q a_k s^k (Durand-Kerner, then Newton polishing).
function polynomial_roots(a::Vector{Float64})
    Q = length(a)
    c = [1.0; a] ./ a[end]
    p(s) = foldr((ck, acc) -> ck + s * acc, c)
    dp(s) = foldr((k, acc) -> k * c[k + 1] + s * acc, 1:Q; init = zero(s))
    z = [(0.4 + 0.9im)^k for k in 1:Q]
    for _ in 1:500
        for i in 1:Q
            d = prod(z[i] - z[j] for j in 1:Q if j != i; init = one(ComplexF64))
            z[i] -= p(z[i]) / d
        end
    end
    for _ in 1:5, i in 1:Q
        z[i] -= p(z[i]) / dp(z[i])
    end
    return z
end

# The roots depend only on (P, Q).
const DENOMINATOR_ROOTS = Dict{Tuple{Int, Int}, Vector{ComplexF64}}()
const DENOMINATOR_ROOTS_LOCK = ReentrantLock()
function denominator_roots(P, Q)
    return lock(DENOMINATOR_ROOTS_LOCK) do
        get!(DENOMINATOR_ROOTS, (P, Q)) do
            w = pade_increment_coefficients(Val(P), Val(Q))[2]
            polynomial_roots([Float64(w[k + 1]) / factorial(k) for k in 1:Q])
        end
    end
end

# An empty object preserves the cache type for :assembled.
function factored_newton(alg, J)
    P, Q = _unwrap_val(alg.order), _unwrap_val(alg.order_q)
    T = eltype(J)
    (is_taylor_pade(alg) && Q >= 3 && T <: AbstractFloat && J isa Matrix) || return nothing
    enabled = alg.newton_matrix in (:factored, :auto)
    return build_factored_newton(
        P, Q, T, enabled ? size(J, 1) : 0;
        auto = alg.newton_matrix === :auto, enabled
    )
end

function build_factored_newton(P, Q, ::Type{T}, n; auto = false, enabled = true) where {T}
    roots = denominator_roots(P, Q)
    real_roots = T[real(r) for r in roots if abs(imag(r)) <= 1.0e-10 * abs(r)]
    complex_roots = Complex{T}[r for r in roots if imag(r) > 1.0e-10 * abs(r)]
    length(real_roots) + 2length(complex_roots) == Q ||
        error("Padé denominator ($P, $Q): unexpected roots $roots")
    return FactoredNewton{T}(
        real_roots, complex_roots,
        [zeros(T, n, n) for _ in real_roots], [zeros(Int, n) for _ in real_roots],
        [zeros(Complex{T}, n, n) for _ in complex_roots], [zeros(Int, n) for _ in complex_roots],
        zeros(Complex{T}, n), zeros(T, n), Ref(T(NaN)),
        enabled, auto, Ref(enabled && !auto), Ref(T(Inf)), Ref(zero(T)), zeros(T, n), zeros(T, n),
        zeros(Complex{T}, n), zeros(Complex{T}, n), Ref(zero(T)), Ref((1, 0, 0))
    )
end

# Estimate of ‖A⁻¹‖₁ from LU factors in LAPACK layout (Hager's method, as in LAPACK's gecon).
# The row permutation of the factorization does not change ‖A⁻¹‖₁ = ‖U⁻¹L⁻¹‖₁, so it is not needed.
function inverse_norm_estimate!(LU, x, y)
    n = size(LU, 1)
    fill!(x, one(eltype(x)) / n)
    estimate = zero(real(eltype(x)))
    for k in 1:5
        copyto!(y, x)
        triangular_solve!(LU, y)
        new_estimate = sum(abs, y)
        k > 1 && new_estimate <= estimate && break
        estimate = new_estimate
        for i in 1:n
            y[i] = iszero(y[i]) ? one(eltype(y)) : y[i] / abs(y[i])
        end
        adjoint_triangular_solve!(LU, y)
        j = argmax(i -> abs(real(y[i])), 1:n)
        abs(real(y[j])) <= real(sum(conj(y[i]) * x[i] for i in 1:n)) && break
        fill!(x, zero(eltype(x)))
        x[j] = one(eltype(x))
    end
    return estimate
end

# y ← U⁻¹ L⁻¹ y (L unit lower, U upper, stored together).
function triangular_solve!(LU, y)
    n = size(LU, 1)
    @inbounds for j in 1:n
        yj = y[j]
        for i in (j + 1):n
            y[i] -= LU[i, j] * yj
        end
    end
    @inbounds for j in n:-1:1
        y[j] /= LU[j, j]
        yj = y[j]
        for i in 1:(j - 1)
            y[i] -= LU[i, j] * yj
        end
    end
    return y
end

# y ← (U⁻¹ L⁻¹)ᴴ y = L⁻ᴴ U⁻ᴴ y.
function adjoint_triangular_solve!(LU, y)
    n = size(LU, 1)
    @inbounds for i in 1:n
        s = y[i]
        for k in 1:(i - 1)
            s -= conj(LU[k, i]) * y[k]
        end
        y[i] = s / conj(LU[i, i])
    end
    @inbounds for i in n:-1:1
        s = y[i]
        for k in (i + 1):n
            s -= conj(LU[k, i]) * y[k]
        end
        y[i] = s
    end
    return y
end

# eps·cond₁(A) from the LU factors of A and ‖A‖₁.
function rounding_floor(LU, anorm, x, y)
    T = real(eltype(LU))
    any(i -> iszero(LU[i, i]), axes(LU, 1)) && return T(Inf)
    return eps(T) * anorm * inverse_norm_estimate!(LU, x, y)
end

# ‖A‖₁, with |re| + |im| ≥ |z| for complex entries (an upper bound within √2, without hypot).
function column_norm(A)
    m = zero(real(eltype(A)))
    @inbounds for j in axes(A, 2)
        s = zero(m)
        for i in axes(A, 1)
            s += abs(real(A[i, j])) + abs(imag(A[i, j]))
        end
        m = max(m, s)
    end
    return m
end

# Contraction rate the step size control counts as fast.
fast_rate(alg) = (c = fast_convergence_cutoff(alg); c / (1 + c))

# Estimate the rounding limit eps·‖M⁻¹(|M||dw|)‖∞/‖dw‖∞ of the assembled solve.
# Switch if it reaches the fast-convergence threshold and the last factored rate.
function check_rounding_limit!(integrator, cache, alg, adaptive, dw)
    F = cache.factored
    (adaptive && F !== nothing && F.enabled && F.auto && cache.lu_valid) || return false
    dnorm = maximum(abs, dw)
    (iszero(dnorm) || !isfinite(dnorm)) && return false
    M, g = cache.Mkept, F.ry
    fill!(g, zero(eltype(g)))
    @inbounds for j in axes(M, 2)
        a = abs(dw[j])
        for i in axes(M, 1)
            g[i] += abs(M[i, j]) * a
        end
    end
    linsolve = cache.linsolve
    u = linsolve.u                     # the correction buffer (dw) stays the solution buffer
    linres = try
        checked_dolinsolve(integrator, linsolve; b = g, linu = F.rx)
    finally
        linsolve.u = u
    end
    integrator.stats.nsolve += 1
    floor = eps(eltype(F.correction)) * maximum(abs, linres.u) / dnorm
    floor >= max(fast_rate(alg), F.rate[]) || return false
    F.active[] = true
    F.distance[] = convert(eltype(F.correction), NaN)
    F.assembled_cap[] = cache.dt_cap
    F.assembled_probe[] = (cache.cap_probe, cache.cap_streak, Int(cache.cap_raised))
    cache.dt_cap = zero(cache.dt_cap)
    cache.cap_probe, cache.cap_streak, cache.cap_raised = 1, 0, false
    return true
end

# Restore the assembled solve and its dt cap when the condition estimate
# or the measured factored rate favors switching.
function return_to_assembled!(cache, F::FactoredNewton, alg, θ)
    F.rate[] = θ
    (F.auto && F.active[] && (F.rounding[] < fast_rate(alg) || θ >= fast_rate(alg))) || return nothing
    F.active[] = false
    cache.dt_cap = F.assembled_cap[]
    probe, streak, raised = F.assembled_probe[]
    cache.cap_probe, cache.cap_streak, cache.cap_raised = probe, streak, raised == 1
    return nothing
end

function reset_factored!(F::FactoredNewton)
    F.distance[] = convert(eltype(F.correction), NaN)
    F.active[] = F.enabled && !F.auto
    F.rounding[] = convert(eltype(F.correction), Inf)
    F.rate[] = zero(eltype(F.correction))
    F.assembled_cap[] = zero(eltype(F.correction))
    F.assembled_probe[] = (1, 0, 0)
    return nothing
end

# In-place LU with partial pivoting (LAPACK layout); false if singular.
function lu_factor!(A, piv)
    n = size(A, 1)
    @inbounds for k in 1:n
        p, m = k, abs(A[k, k])
        for i in (k + 1):n
            a = abs(A[i, k])
            a > m && ((p, m) = (i, a))
        end
        piv[k] = p
        (iszero(m) || !isfinite(m)) && return false
        if p != k
            for j in 1:n
                A[k, j], A[p, j] = A[p, j], A[k, j]
            end
        end
        inv_pivot = inv(A[k, k])
        for i in (k + 1):n
            A[i, k] *= inv_pivot
        end
        for j in (k + 1):n
            akj = A[k, j]
            iszero(akj) && continue
            for i in (k + 1):n
                A[i, j] -= A[i, k] * akj
            end
        end
    end
    return true
end

function lu_solve!(A, piv, b)
    n = size(A, 1)
    @inbounds for k in 1:n
        p = piv[k]
        p != k && ((b[k], b[p]) = (b[p], b[k]))
    end
    @inbounds for j in 1:n
        bj = b[j]
        iszero(bj) && continue
        for i in (j + 1):n
            b[i] -= A[i, j] * bj
        end
    end
    @inbounds for j in n:-1:1
        b[j] /= A[j, j]
        bj = b[j]
        for i in 1:(j - 1)
            b[i] -= A[i, j] * bj
        end
    end
    return b
end

# Factor I - (d/r_i) J for every root and record eps·Π cond(factor); false if a factor is singular.
function factor_denominator!(F::FactoredNewton, Jf, d)
    n = size(Jf, 1)
    T = eltype(F.correction)
    rounding = eps(T)
    for (k, r) in enumerate(F.real_roots)
        A = F.real_factors[k]
        s = d / r
        @inbounds for j in 1:n, i in 1:n
            A[i, j] = -s * Jf[i, j]
        end
        @inbounds for i in 1:n
            A[i, i] += one(eltype(A))
        end
        anorm = F.auto ? column_norm(A) : zero(T)
        lu_factor!(A, F.real_pivots[k]) || return false
        F.auto && (rounding *= rounding_floor(A, anorm, F.rx, F.ry) / eps(T))
    end
    for (k, r) in enumerate(F.complex_roots)
        A = F.complex_factors[k]
        s = d / r
        @inbounds for j in 1:n, i in 1:n
            A[i, j] = -s * Jf[i, j]
        end
        @inbounds for i in 1:n
            A[i, i] += one(eltype(A))
        end
        anorm = F.auto ? column_norm(A) : zero(T)
        lu_factor!(A, F.complex_pivots[k]) || return false
        # The conjugate factor has the same condition number.
        F.auto && (rounding *= (rounding_floor(A, anorm, F.cx, F.cy) / eps(T))^2)
    end
    F.rounding[] = rounding
    return true
end

# b ← D(dJ)^{-1} b, one factor after the other.
function solve_denominator!(F::FactoredNewton, b)
    for k in eachindex(F.real_roots)
        lu_solve!(F.real_factors[k], F.real_pivots[k], b)
    end
    for k in eachindex(F.complex_roots)
        A, piv, y = F.complex_factors[k], F.complex_pivots[k], F.cbuf
        y .= b
        lu_solve!(A, piv, y)          # (I - (d/r) J)^{-1}
        conj!(y)
        lu_solve!(A, piv, y)          # conj((I - (d/r̄) J)^{-1} conj(y))
        conj!(y)
        b .= real.(y)
    end
    return b
end
nfactors(F::FactoredNewton) = length(F.real_roots) + length(F.complex_roots)
nsolves(F::FactoredNewton) = length(F.real_roots) + 2length(F.complex_roots)

# The stage solve with the factored matrix. Same stopping rule, validation and step size
# feedback as the assembled iteration (see implicit_taylor_stage!).
function factored_stage!(integrator, cache, F::FactoredNewton, t, dt, uprev, u, tc, distance, weights)
    (; κ, tmp, tmp_embedded, atmp, uintermediate, rhs, utilde) = cache
    alg = unwrap_alg(integrator, true)
    (; internalnorm, abstol, reltol) = integrator.opts
    stats = integrator.stats
    Jf = cache.Mtmp            # ∂f/∂u at the last refresh (Jcheck aliases it)
    dz = F.correction
    taylor_initial_guess!(integrator, cache, alg, t, dt, uprev, true)
    @. uintermediate = uprev + utilde
    ndw = ndwprev = θ = η = one(κ)
    θ_slow = convert(typeof(κ), 1 // 2)
    κ_second = κ * convert(typeof(κ), 1 // 10)
    tiny = κ * convert(typeof(κ), 1 // 1000)
    # Factors are kept while the distance changes by less than the matrix cutoff.
    refresh = !(abs(distance - F.distance[]) <= abs(distance) * matrix_cutoff(alg))
    converged = false
    measured = false
    fresh_iter = 1
    restarted = false
    validation_failed = false
    iter = 0
    while iter < alg.maxiters
        iter += 1
        stats.nnonliniter += 1
        if refresh
            cache.polynomial_check(tmp, tmp_embedded, cache.Jcheck, uintermediate, tc, distance)
            stats.njacs += 1
            cache.ncheckjac += 1
            stats.nw += nfactors(F)
            F.distance[] = convert(eltype(F.correction), NaN)
            factor_denominator!(F, Jf, distance) || break
            F.distance[] = distance
            fresh_iter = iter
        else
            cache.polynomial_iter(tmp, tmp_embedded, uintermediate, tc, distance)
        end
        dz .= _vec(utilde) .+ _vec(tmp) .- _vec(rhs)
        solve_denominator!(F, dz)
        stats.nsolve += nsolves(F)
        iter > 1 && (ndwprev = ndw)
        correction = shaped_like(dz, utilde)
        @.. broadcast = false atmp = calculate_residuals(
            abs(correction), uprev, uintermediate, abstol, reltol, internalnorm, t
        )
        ndw = internalnorm(atmp, t)
        isfinite(ndw) || break
        @. utilde = utilde - correction
        @. uintermediate = uprev + utilde
        measured = iter >= fresh_iter + 1
        if measured
            θ = ndw / ndwprev
            η = θ / (1 - θ)
        end
        if ndw <= tiny || (iter == fresh_iter + 1 && ndw < κ_second) ||
                (measured && θ < one(θ) && η * ndw < κ)
            defect, _, paired = check_stage_defect!(integrator, cache, t, tc, distance, uprev, weights)
            if isfinite(defect) && defect <= κ
                cache.embedded_filtered = paired
                converged = true
                break
            end
            cache.nvalidation_fail += 1
            (validation_failed || !isfinite(defect)) && break
            validation_failed = true
            cache.recovery_since_accept = true
            refresh = true
            continue
        end
        if measured && θ > 2
            # Diverging: start over from the predictor with factors formed there, once.
            restarted && break
            restarted = true
            cache.nrestart += 1
            cache.recovery_since_accept = true
            taylor_initial_guess!(integrator, cache, alg, t, dt, uprev, true)
            @. uintermediate = uprev + utilde
            refresh = true
            continue
        end
        refresh = measured && θ > θ_slow
    end
    cache.iter = iter
    if !converged
        F.distance[] = convert(eltype(F.correction), NaN)
        cache.jac_valid = false
        cache.newton_failed = true
        record_newton_failure!(cache, dt)
        integrator.force_stepfail = true
        stats.nnonlinconvfail += 1
        cache.status = Divergence
        return false
    end
    cache.status = Convergence
    record_newton_convergence!(cache, measured, η, alg)
    return_to_assembled!(cache, F, alg, measured ? convert(eltype(F.correction), θ) : zero(eltype(F.correction)))
    u .= uintermediate
    return true
end
