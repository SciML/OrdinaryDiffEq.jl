import OrdinaryDiffEqTaylorSeries: build_polynomial_bundles, build_stage_bundles

# The polynomial kernels evaluate weighted sums of a Taylor jet. polynomial_prev
# forms the explicit side and predictor. polynomial_iter forms the implicit side
# and error polynomial. polynomial_jac and polynomial_check also write ∂c_k/∂u
# for k = 1:Q and k = 1, respectively (degree P for μ-Taylor).
# stage_values writes the jet of f along the free stages. stage_jac also writes f_u.
@cache mutable struct ImplicitTaylorCache{
        T, uType, tType, stageTimeType, jacobianType, jacobianStackType, rowBufferType, uIntermediateType,
        rateType, uNoUnitsType, F1, F2, Tol, FN, stageType, stageJacType, SF,
    } <:
    OrdinaryDiffEqMutableCache
    μ::T
    t::tType               # start time of the current step
    polynomial_prev::FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, uIntermediateType, uIntermediateType, stageTimeType, stageTimeType}}
    polynomial_iter::FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, uIntermediateType, stageTimeType, stageTimeType}}
    polynomial_jac::FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, jacobianStackType, uIntermediateType, stageTimeType, stageTimeType}}
    polynomial_check::FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, jacobianStackType, uIntermediateType, stageTimeType, stageTimeType}}
    stage_values::FunctionWrapper{Nothing, Tuple{stageType, uIntermediateType, stageType, stageTimeType, stageTimeType}}
    stage_jac::FunctionWrapper{Nothing, Tuple{stageType, stageJacType, uIntermediateType, stageType, stageTimeType, stageTimeType}}
    u::uType
    uprev::uType
    uprev2::uType
    # State and time two steps back (extrapolation, third divided differences).
    uprev3::uType
    uhist::uType           # uprev2 of the previous step, becomes uprev3 at the next
    tprev2::tType
    t_hist::tType
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
    # One factor W of the error filter, its scaled Jacobian JW, and its solver.
    # For degree one the solver and matrix alias those of the Newton solve.
    W::jacobianType
    JW::jacobianType
    linsolve_filter::F2
    filter_valid::Bool
    embedded_filtered::Bool
    # The solver returns its solution in its own buffer instead of the cache's `u`.
    detached_solution::Bool
    detached_filter::Bool
    direct_lu::Bool            # Newton matrix and filter: dense_lu!/dense_solve! in place
    pivJ::Vector{Int}
    pivW::Vector{Int}
    κ::Tol
    ηold::Tol
    ηleak::Tol         # last measured stiff rate of the leak test (leak_stop_rate) ...
    ηleak_h::tType     # ... and the step size it was measured at
    status::NLStatus
    iter::Int          # total corrections in this stage solve, including a restart
    ncoefjac::Int      # full coefficient-Jacobian evaluations
    ncheckjac::Int     # endpoint validation Jacobian evaluations (including degree 1)
    nendpoint_refresh::Int
    nrestart::Int
    nvalidation_fail::Int
    nfilter_recovery::Int # bounded filtered-direction recovery attempts
    nleak::Int         # stage candidates rejected by the leak test (derivative stages)
    reduced_check_h::stageTimeType # last distance where only the lifted check passed
    jac_valid::Bool             # whether Jk has been formed
    constant_jacobian::Bool     # symbolic proof that every Jk is independent of (u, t)
    lu_valid::Bool              # whether J holds a factorization (of Mkept)
    newton_failed::Bool # a Newton iteration failed since the last accepted step
    newton_limit::Bool  # slow convergence or recovery: do not grow dt
    newton_rate_measured::Bool # whether the final attempt measured a contraction rate
    newton_slow::Bool   # measured slow convergence, independent of recovery
    recovery_since_accept::Bool # nonlinear recovery; retained through error retries until acceptance
    dt_cap::tType       # |dt| above which the Newton iteration last failed (0: none)
    # Fast steps required at the step-size cap before lifting it and probing again.
    cap_probe::Int
    cap_streak::Int
    cap_raised::Bool
    factored::FN        # FactoredNewton (newton_matrix = :auto or :factored), otherwise nothing
    # V and Fst are n×Q, and Jst holds the single n×n endpoint Jacobian.
    stage_system::Bool
    # Automatic switching probes the reduced form after fast derivative-stage steps.
    stage_auto::Bool
    use_stages::Bool
    stage_trial::Bool
    stage_probe::Int
    stage_streak::Int
    stage_recovery::Bool
    stage_handover::Bool # a same-step rescue; retained through error retries until acceptance
    V::stageType
    V0::stageType
    stages_h::stageTimeType
    stages_valid::Bool
    Fst::stageType
    Jst::stageJacType
    stage_factors::SF              # factors (I + λ_i hJ) of the simplified Newton matrix
    rb::uIntermediateType          # right-hand side (Qn) and solution dv of the stage system
    dv::uIntermediateType
    r0::uIntermediateType
    dz::uIntermediateType
    # Start-of-step Jacobian and scratch space for stage_leak_norm!.
    Jprev::stageJacType
    jprev_valid::Bool
    leak_p::uIntermediateType
    leak_q::uIntermediateType
    leak_w::uIntermediateType
    stage_jac_valid::Bool
    constant_stage_jacobian::Bool
    lu_stages_valid::Bool
end

function alg_cache(
        alg::ImplicitTaylor{P, Q}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {P, Q, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return implicit_taylor_cache(
        alg, u, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, uprev, uprev2, f, t, dt, p
    )
end

function implicit_taylor_cache(
        alg::ImplicitTaylor{P, Q}, u, rate_prototype,
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
            # Order-two methods use this until divided differences have enough states.
            P1, Q1 = embedded_degrees(P, Q)
            tuple_A1, tuple_B1 = pade_increment_coefficients(Val(P1), Val(Q1))
        else
            # Backward Euler: no embedded formula; the second outputs are unused.
            tuple_A1, tuple_B1 = tuple_A, tuple_B
        end
        implicit_specs = (
            (coefficient_sets = (tuple_B, tuple_B1), coefficient_jacobians = 0),
            (coefficient_sets = (tuple_B, tuple_B1), coefficient_jacobians = Q),
            (coefficient_sets = (tuple_B, tuple_B1), coefficient_jacobians = 1),
        )
        if P == 0
            polynomial_prev = nothing
            polynomial_iter, polynomial_jac, polynomial_check =
                build_polynomial_bundles(polynomial_f, Val(true), p, implicit_specs, n)
        else
            polynomial_start = alg.extrapolant === :polynomial && P + Q > 2
            previous_sets = polynomial_start ? (tuple_A, tuple_A1, tuple_B) : (tuple_A, tuple_A1)
            previous, polynomial_iter, polynomial_jac, polynomial_check =
                build_polynomial_bundles(
                polynomial_f, Val(true), p,
                ((coefficient_sets = previous_sets, coefficient_jacobians = 0), implicit_specs...), n
            )
            # Only the polynomial predictor needs B_inc(uprev).
            polynomial_prev = polynomial_start ? previous :
                (a, a1, b0, u, t, h) -> previous(a, a1, u, t, h)
        end
        degree = Q
    else
        # μ-Taylor: z + T_inc(uprev + z, t + μh, -μh) = 0, then
        # u = u_c + T_inc(u_c, t_c, (1 - μ)h); c_P h^P gives the error estimate.
        tuple_T = taylor_increment_coefficients(Val(P))
        tuple_L = taylor_last_coefficients(Val(P))
        polynomial_prev = nothing
        polynomial_iter, polynomial_jac, polynomial_check =
            build_polynomial_bundles(
            polynomial_f, Val(true), p, (
                (coefficient_sets = (tuple_T, tuple_L), coefficient_jacobians = 0),
                (coefficient_sets = (tuple_T, tuple_L), coefficient_jacobians = P),
                (coefficient_sets = (tuple_T, tuple_L), coefficient_jacobians = 1),
            ), n
        )
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
    # Cache row norms for the column-major scans at n > 32.
    row_change = zeros(typeof(real(zero(eltype(J)))), n > 32 ? n : 0)
    Mrow_norm = n > 32 ? similar(row_change) : row_change
    JWrow_norm = n > 32 ? similar(row_change) : row_change
    # Small dense assembled solves use the in-place LU kernels by default.
    factored = factored_newton(alg, J)
    direct_lu = alg.linsolve === nothing && J isa Matrix && n <= DIRECT_LU_MAX &&
        alg.newton_matrix === :assembled
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
    detached_solution = direct_lu ? false : solution_detached(linsolve, J)
    detached_filter = direct_lu ? false :
        degree > 1 ? solution_detached(linsolve_filter, W) : detached_solution
    pivJ = zeros(Int, direct_lu ? n : 0)
    pivW = zeros(Int, direct_lu && degree > 1 ? n : 0)

    # Pairing conjugate Padé roots requires real stage equations. The factored
    # Newton matrix option belongs to the increment form.
    stage_system = is_taylor_pade(alg) && Q >= 2 && alg.stages in (:derivative, :auto) &&
        alg.newton_matrix === :assembled && J isa Matrix && eltype(J) <: Real
    uIntermediateType = typeof(tmp)
    if stage_system
        stage_values, stage_jac = build_stage_bundles(
            polynomial_f, Val(true), p, Q, n;
            jacobian_order = 0
        )
        V, V0, Fst = zeros(eltype(tmp), n, Q), zeros(eltype(tmp), n, Q), zeros(eltype(tmp), n, Q)
        Jst = zeros(eltype(J), n, n, 1)
        stage_factors = build_stage_factors(
            tuple_B, J, alg.linsolve, linear_u, p, t;
            direct = alg.linsolve === nothing && n <= DIRECT_LU_MAX
        )
        rb, dv = zeros(eltype(tmp), Q * n), zeros(eltype(tmp), Q * n)
        r0, dz = similar(tmp), similar(tmp)
        Jprev = zeros(eltype(J), n, n, 1)
        leak_p, leak_q, leak_w = similar(tmp), similar(tmp), similar(tmp)
        constant_stage_jacobian = stage_jac.constant_jacobian
    else
        stage_values = (F, u, V, t, h) -> nothing
        stage_jac = (F, J, u, V, t, h) -> nothing
        V, V0, Fst = zeros(eltype(tmp), n, 0), zeros(eltype(tmp), n, 0), zeros(eltype(tmp), n, 0)
        Jst = zeros(eltype(J), n, n, 0)
        # Empty factors preserve the cache type when the stage system is disabled.
        stage_factors = empty_stage_factors(real(eltype(J)), alg.linsolve, linear_u, p, t)
        rb, dv = tmp_vec, tmp_vec
        r0, dz = tmp, tmp
        Jprev = Jst
        leak_p, leak_q, leak_w = tmp, tmp, tmp
        constant_stage_jacobian = false
    end
    stageType, stageJacType = typeof(V), typeof(Jst)
    valuestype = FunctionWrapper{Nothing, Tuple{stageType, uIntermediateType, stageType, stageTimeType, stageTimeType}}
    stagejactype = FunctionWrapper{Nothing, Tuple{stageType, stageJacType, uIntermediateType, stageType, stageTimeType, stageTimeType}}

    ηold = one(uToltype)
    iter = 0
    status = Convergence
    jacobianType, jacobianStackType = typeof(J), typeof(Jk)
    prevtype = FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, uIntermediateType, uIntermediateType, stageTimeType, stageTimeType}}
    itertype = FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, uIntermediateType, stageTimeType, stageTimeType}}
    jactype = FunctionWrapper{Nothing, Tuple{uIntermediateType, uIntermediateType, jacobianStackType, uIntermediateType, stageTimeType, stageTimeType}}
    wrapped_prev = polynomial_prev === nothing ?
        prevtype((a, a1, b0, u, t, h) -> nothing) : prevtype(polynomial_prev)

    return ImplicitTaylorCache(
        alg.μ, t, wrapped_prev, itertype(polynomial_iter), jactype(polynomial_jac), jactype(polynomial_check),
        valuestype(stage_values), stagejactype(stage_jac),
        u, uprev, uprev2, zero(u), zero(u), t, t, J, Jk, Mkept, Mtmp, Jcheck,
        Mrow_norm, JWrow_norm, row_change, utilde, uintermediate, rhs, rhs_embedded, predictor,
        tmp, tmp_embedded, stage_error, atmp, fsalfirst, fsallast, linsolve, W, JW, linsolve_filter, false, false,
        detached_solution, detached_filter, direct_lu, pivJ, pivW,
        κ, ηold, ηold, zero(t), status, iter, 0, 0, 0, 0, 0, 0, 0, zero(stageTimeType),
        false, polynomial_jac.constant_jacobian, false, false, false, false, false, false, zero(t), 1, 0, false, factored,
        stage_system, stage_system && alg.stages === :auto, stage_system && alg.stages === :derivative,
        false, 1, 0, false, false,
        V, V0, zero(stageTimeType), false, Fst, Jst, stage_factors,
        rb, dv, r0, dz, Jprev, false, leak_p, leak_q, leak_w,
        false, constant_stage_jacobian, false
    )
end

struct TaylorLinearSolveFailure <: Exception
    message::String
end

Base.showerror(io::IO, err::TaylorLinearSolveFailure) = print(io, err.message)

@inline function check_finite_solution(u)
    all(isfinite, u) ||
        throw(TaylorLinearSolveFailure("ImplicitTaylor linear solve returned a nonfinite solution"))
    return u
end

@inline function check_linear_solution(result)
    SciMLBase.successful_retcode(result.retcode) ||
        throw(TaylorLinearSolveFailure("ImplicitTaylor linear solve failed: $(result.retcode)"))
    check_finite_solution(result.u)
    return result
end

@inline checked_dolinsolve(args...; kwargs...) = check_linear_solution(dolinsolve(args...; kwargs...))

# Probe whether the linear solver writes into its supplied solution buffer.
function solution_detached(linsolve, A)
    fill!(A, zero(eltype(A)))
    @inbounds for i in 1:size(A, 1)
        A[i, i] = one(eltype(A))
    end
    LinearSolve.reinit!(linsolve; A)
    fill!(linsolve.b, one(eltype(linsolve.b)))
    result = check_linear_solution(LinearSolve.solve!(linsolve))
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
    cache = implicit_taylor_cache(
        alg, state, rates, uEltypeNoUnits,
        uBottomEltypeNoUnits, copy(state), copy(state), SciMLBase.ODEFunction{true}(f!), t, dt, p; linear_u = u
    )
    return ImplicitTaylorConstantCache(cache)
end

@inline function SciMLBase.get_tmp_cache(
        integrator, alg::Union{ImplicitTaylor},
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp, cache.atmp)
end
