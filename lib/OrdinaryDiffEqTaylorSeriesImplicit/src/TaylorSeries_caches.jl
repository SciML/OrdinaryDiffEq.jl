import OrdinaryDiffEqTaylorSeries: build_polynomial_bundle

# All quantities are weighted sums of the same Taylor jet coefficients,
# evaluated by four generated kernels:
#   polynomial_prev(a, a1, b0, uprev, t, h)   explicit side at the old state
#   polynomial_iter(b, b1, u, t, h)           implicit side at a Newton iterate
#   polynomial_jac(b, b1, Jk, u, t, h)        same, plus ∂c_k/∂u, k = 1:Q
#   polynomial_check(b, b1, J1, u, t, h)      same, plus only ∂c_1/∂u
# b1 is the embedded implicit polynomial (Taylor-Padé) or the highest Taylor
# term c_P h^P (μ-Taylor). The ∂c_k/∂u do not depend on h, so the Newton matrix
# M = I + Σ_k w_k h^k ∂c_k/∂u can be assembled from them for any h.
@cache mutable struct ImplicitTaylorCache{
        T, uType, tType, stageTimeType, jacobianType, jacobianStackType, rowBufferType, uIntermediateType,
        rateType, uNoUnitsType, F1, F2, Tol, StepLimiter,
    } <:
    OrdinaryDiffEqMutableCache
    μ::T
    t::tType               # start time of the current step
    polynomial_prev::FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, uIntermediateType, uIntermediateType, stageTimeType, stageTimeType}}
    polynomial_iter::FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, uIntermediateType, stageTimeType, stageTimeType}}
    polynomial_jac::FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, jacobianStackType, uIntermediateType, stageTimeType, stageTimeType}}
    polynomial_check::FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, jacobianStackType, uIntermediateType, stageTimeType, stageTimeType}}
    u::uType
    uprev::uType
    uprev2::uType
    # State and time two steps back (extrapolation, third divided differences).
    uprev3::uType
    uhist::uType           # uprev2 of the previous step, becomes uprev3 at the next
    tprev2::tType
    thist::tType
    J::jacobianType        # Newton matrix M (holds its factorization after the solve)
    Jk::jacobianStackType  # ∂c_k/∂u, k = 1:Q, from the last Jacobian evaluation
    Mkept::jacobianType             # matrix whose factorization J holds
    Mtmp::jacobianType              # assembly scratch; holds endpoint f' after validation
    Jcheck::jacobianStackType       # singleton third-dimension view of Mtmp
    Mrow_norm::rowBufferType
    JWrow_norm::rowBufferType
    row_change::rowBufferType
    utilde::uIntermediateType        # Newton increment z (stage value - uprev)
    uintermediate::uIntermediateType # uprev + z
    rhs::uIntermediateType           # A_inc(uprev): explicit right-hand side
    rhs_embedded::uIntermediateType  # A1_inc(uprev): embedded right-hand side
    predictor::uIntermediateType     # B_inc(uprev): for the polynomial predictor
    tmp::uIntermediateType           # B_inc at the iterate, then the residual
    tmp_embedded::uIntermediateType  # B1_inc at the iterate / highest Taylor term
    stage_error::uIntermediateType   # filtered main residual for convergence validation
    atmp::uNoUnitsType
    fsalfirst::rateType
    fsallast::rateType
    linsolve::F1
    # Error filter (I + (w_1/Q) d ∂c_1/∂u)^Q: the factor W, the scaled Jacobian
    # JW it was formed from, and its solver. For Q = 1 these alias M.
    W::jacobianType
    JW::jacobianType
    linsolve_filter::F2
    filter_valid::Bool
    embedded_filtered::Bool
    # The solver returns its solution in its own buffer instead of the cache's `u`.
    detached_solution::Bool
    detached_filter::Bool
    κ::Tol
    ηold::Tol
    status::NLStatus
    iter::Int          # total corrections in this stage solve, including a restart
    ncoefjac::Int      # full coefficient-Jacobian evaluations
    ncheckjac::Int     # endpoint validation Jacobian evaluations (including degree 1)
    nendpoint_refresh::Int
    nrestart::Int
    nvalidation_fail::Int
    nfilter_recovery::Int # bounded filtered-direction recovery attempts
    jac_valid::Bool             # whether Jk has been formed
    lu_valid::Bool              # whether J holds a factorization (of Mkept)
    newton_failed::Bool # a Newton iteration failed since the last accepted step
    newton_limit::Bool  # the last Newton iteration converged slowly: do not grow dt
    recovery_since_accept::Bool # retained across rejected attempts, consumed on acceptance
    dt_cap::tType       # |dt| above which the Newton iteration last failed (0: none)
    # Cap probing (see newton_feedback!): fast steps needed at the cap before it
    # is raised, fast steps so far, and whether it was just raised.
    cap_probe::Int
    cap_streak::Int
    cap_raised::Bool
    step_limiter!::StepLimiter
end

function alg_cache(
        alg::ImplicitTaylor{P, Q}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {P, Q, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return implicit_taylor_cache(alg, u, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, uprev, uprev2, f, t, dt, p)
end

function implicit_taylor_cache(alg::ImplicitTaylor{P, Q}, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, uprev, uprev2, f, t, dt, p;
        linear_u = u
    ) where {P, Q, uEltypeNoUnits, uBottomEltypeNoUnits}
    # Unwrap the ODEFunction once for all generated kernels.
    polynomial_f = SciMLBase.unwrapped_f(f)
    n = length(u)
    if is_taylor_pade(alg)
        # Increment form: z + B_inc(uprev + z, t + h, h) = A_inc(uprev, t, h).
        tuple_A, tuple_B = pade_increment_coefficients(Val(P), Val(Q))
        if uses_embedded(alg) || uses_divided3(alg)
            # Embedded formula of one order less (see alg_utils.jl); second-order
            # methods use it only until the divided differences have enough states.
            P1, Q1 = embedded_degrees(P, Q)
            tuple_A1, tuple_B1 = pade_increment_coefficients(Val(P1), Val(Q1))
        else
            # Backward Euler: no embedded formula; the second outputs are unused.
            tuple_A1, tuple_B1 = tuple_A, tuple_B
        end
        polynomial_prev = P > 0 ?
            build_polynomial_bundle(polynomial_f, Val(true), p, (tuple_A, tuple_A1, tuple_B), n) :
            nothing
        polynomial_iter = build_polynomial_bundle(polynomial_f, Val(true), p, (tuple_B, tuple_B1), n)
        polynomial_jac = build_polynomial_bundle(polynomial_f, Val(true), p, (tuple_B, tuple_B1), n;
            coefficient_jacobians = Q)
        polynomial_check = build_polynomial_bundle(polynomial_f, Val(true), p,
            (tuple_B, tuple_B1), n; coefficient_jacobians = 1)
        degree = Q
    else
        # μ-Taylor: z + T_inc(uprev + z, t + μh, -μh) = 0, then
        # u = u_c + T_inc(u_c, t_c, (1 - μ)h); c_P h^P gives the error estimate.
        tuple_T = taylor_increment_coefficients(Val(P))
        tuple_L = taylor_last_coefficients(Val(P))
        polynomial_prev = nothing
        polynomial_iter = build_polynomial_bundle(polynomial_f, Val(true), p, (tuple_T, tuple_L), n)
        polynomial_jac = build_polynomial_bundle(polynomial_f, Val(true), p, (tuple_T, tuple_L), n;
            coefficient_jacobians = P)
        polynomial_check = build_polynomial_bundle(polynomial_f, Val(true), p,
            (tuple_T, tuple_L), n; coefficient_jacobians = 1)
        degree = P
    end
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    # Residual norms are real even when the internal expansion center is complex.
    atmp = similar(u, typeof(real(zero(uEltypeNoUnits))))
    recursivefill!(atmp, false)
    uToltype = constvalue(typeof(real(zero(uBottomEltypeNoUnits))))
    promoted_eltype = is_mu_taylor(alg) ? promote_type(eltype(u), typeof(alg.μ)) : eltype(u)
    stageTimeType = is_mu_taylor(alg) ? promote_type(typeof(t), typeof(alg.μ)) : typeof(t)
    tmp = similar(u, promoted_eltype)
    tmp_embedded = similar(tmp)
    stage_error = similar(tmp)
    rhs, rhs_embedded, predictor = similar(tmp), similar(tmp), similar(tmp)
    uintermediate = similar(tmp)
    utilde = similar(tmp)
    tmp_vec = _vec(tmp)
    linsolve_u0 = copy(tmp_vec)
    κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1 // 100)
    J = ArrayInterface.zeromatrix(tmp_vec)
    Mkept, Mtmp = ArrayInterface.zeromatrix(tmp_vec), ArrayInterface.zeromatrix(tmp_vec)
    Jcheck = reshape(Mtmp, n, n, 1)
    Jk = zeros(eltype(J), n, n, degree)
    # n > 32: column-wise row scans with cached row norms; otherwise an early-exit
    # row scan (empty buffers keep the cache type independent of n).
    row_change = zeros(typeof(real(zero(eltype(J)))), n > 32 ? n : 0)
    Mrow_norm = n > 32 ? similar(row_change) : row_change
    JWrow_norm = n > 32 ? similar(row_change) : row_change
    linprob = LinearProblem(J, tmp_vec, (nothing, linear_u, p, t); u0 = linsolve_u0)
    linsolve = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        assumptions = LinearSolve.OperatorAssumptions(true)
    )
    if degree > 1
        W = ArrayInterface.zeromatrix(tmp_vec)
        JW = ArrayInterface.zeromatrix(tmp_vec)
        # Default solvers filter the main and embedded defects as two columns.
        filter_rhs = alg.linsolve === nothing && uses_embedded(alg) ?
            zeros(eltype(tmp), length(tmp), 2) : copy(tmp_vec)
        filter_prob = LinearProblem(W, filter_rhs, (nothing, linear_u, p, t); u0 = copy(filter_rhs))
        filter_alg = alg.linsolve
        if filter_alg === nothing
            # Choose by the state dimension, as for a vector right-hand side.
            filter_alg = LinearSolve.defaultalg(W, tmp_vec, LinearSolve.OperatorAssumptions(true))
            # Small filters: generic LU avoids the call overhead of Accelerate.
            if 10 < n <= 24 && eltype(W) <: Union{Float32, Float64} &&
                    filter_alg.alg === LinearSolve.DefaultAlgorithmChoice.AppleAccelerateLUFactorization
                filter_alg = LinearSolve.DefaultLinearSolver(LinearSolve.DefaultAlgorithmChoice.GenericLUFactorization)
            end
        end
        linsolve_filter = init(
            filter_prob, filter_alg, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
            assumptions = LinearSolve.OperatorAssumptions(true)
        )
    else
        W, JW, linsolve_filter = J, J, linsolve
    end
    detached_solution = solution_detached(linsolve, J)
    detached_filter = degree > 1 ? solution_detached(linsolve_filter, W) : detached_solution

    ηold = one(uToltype)
    iter = 0
    status = Convergence
    uIntermediateType, jacobianType, jacobianStackType = typeof(tmp), typeof(J), typeof(Jk)
    prevtype = FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, uIntermediateType, uIntermediateType, stageTimeType, stageTimeType}}
    itertype = FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, uIntermediateType, stageTimeType, stageTimeType}}
    jactype = FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, jacobianStackType, uIntermediateType, stageTimeType, stageTimeType}}
    wrapped_prev = polynomial_prev === nothing ?
        prevtype((a, a1, b0, u, t, h) -> nothing) : prevtype(polynomial_prev)

    return ImplicitTaylorCache(
        alg.μ, t, wrapped_prev, itertype(polynomial_iter), jactype(polynomial_jac), jactype(polynomial_check),
        u, uprev, uprev2, zero(u), zero(u), t, t, J, Jk, Mkept, Mtmp, Jcheck,
        Mrow_norm, JWrow_norm, row_change, utilde, uintermediate, rhs, rhs_embedded, predictor,
        tmp, tmp_embedded, stage_error, atmp, fsalfirst, fsallast, linsolve, W, JW, linsolve_filter, false, false,
        detached_solution, detached_filter,
        κ, ηold, status, iter, 0, 0, 0, 0, 0, 0,
        false, false, false, false, false, zero(t), 1, 0, false, alg.step_limiter!
    )
end

# Solve I x = 1 once to see whether `solve!` writes into the cache's `u` or
# returns its own buffer. A is scratch here and refactorized before use.
function solution_detached(linsolve, A)
    fill!(A, zero(eltype(A)))
    @inbounds for i in 1:size(A, 1)
        A[i, i] = one(eltype(A))
    end
    LinearSolve.reinit!(linsolve; A)
    fill!(linsolve.b, one(eltype(linsolve.b)))
    result = LinearSolve.solve!(linsolve)
    return result.u !== linsolve.u
end

function get_fsalfirstlast(cache::ImplicitTaylorCache, u)
    return (cache.fsalfirst, cache.fsallast)
end

struct ImplicitTaylorConstantCache{C} <: OrdinaryDiffEqConstantCache
    stage_cache::C
end

# Out-of-place problems are solved through the in-place vector stage equations.
taylor_state(u::Number, v) = only(v)
taylor_state(u::AbstractArray, v) = reshape(v, size(u))
taylor_result(u::Number, v) = only(v)
taylor_result(u::AbstractArray, v) = ArrayInterface.restructure(u, copy(v))

function alg_cache(
        alg::ImplicitTaylor, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2,
        f, t, dt, reltol, p, calck, ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    state = u isa Number ? [u] : collect(vec(u))
    rates = rate_prototype isa Number ? [rate_prototype] : collect(vec(rate_prototype))
    original_f = SciMLBase.unwrapped_f(f)
    f! = function (du, v, p, t)
        value = original_f(taylor_state(u, v), p, t)
        if value isa Number
            du[1] = value
        else
            copyto!(du, vec(value))
        end
        return nothing
    end
    cache = implicit_taylor_cache(alg, state, rates, uEltypeNoUnits,
        uBottomEltypeNoUnits, copy(state), copy(state), SciMLBase.ODEFunction{true}(f!), t, dt, p; linear_u = u)
    return ImplicitTaylorConstantCache(cache)
end

@inline function SciMLBase.get_tmp_cache(
        integrator, alg::Union{ImplicitTaylor},
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp, cache.atmp)
end
