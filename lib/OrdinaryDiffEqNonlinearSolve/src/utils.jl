get_status(nlsolver::AbstractNLSolver) = nlsolver.status
get_new_W_γdt_cutoff(nlsolver::AbstractNLSolver) = nlsolver.cache.new_W_γdt_cutoff
# handle FIRK
get_new_W_γdt_cutoff(alg::NewtonAlgorithm) = alg.new_W_γdt_cutoff

nlsolvefail(nlsolver::AbstractNLSolver) = nlsolvefail(get_status(nlsolver))
nlsolvefail(status::NLStatus) = Int8(status) <= 0

relax(nlsolver::AbstractNLSolver) = relax(nlsolver.alg)
relax(alg::NLNewton) = alg.relax
relax(_) = 0 // 1

isnewton(nlsolver::AbstractNLSolver) = isnewton(nlsolver.cache)
isnewton(::AbstractNLSolverCache) = false
isnewton(::Union{NLNewtonCache, NLNewtonConstantCache}) = true

is_always_new(nlsolver::AbstractNLSolver) = is_always_new(nlsolver.alg)
check_div(nlsolver::AbstractNLSolver) = check_div(nlsolver.alg)
check_div(alg) = isdefined(alg, :check_div) ? alg.check_div : true

isJcurrent(nlsolver::AbstractNLSolver, integrator) = integrator.t == nlsolver.cache.J_t
isfirstcall(nlsolver::AbstractNLSolver) = nlsolver.cache.firstcall
isfirststage(nlsolver::AbstractNLSolver) = nlsolver.cache.firststage
setfirststage!(nlsolver::AbstractNLSolver, val::Bool) = setfirststage!(nlsolver.cache, val)
function setfirststage!(nlcache::Union{NLNewtonCache, NLNewtonConstantCache}, val::Bool)
    (nlcache.firststage = val)
end
setfirststage!(::Any, val::Bool) = nothing
markfirststage!(nlsolver::AbstractNLSolver) = setfirststage!(nlsolver, true)

getnfails(_) = 0
getnfails(nlsolver::AbstractNLSolver) = nlsolver.nfails

set_new_W!(nlsolver::AbstractNLSolver, val::Bool)::Bool = set_new_W!(nlsolver.cache, val)
set_new_W!(nlcache::Union{NLNewtonCache, NLNewtonConstantCache}, val::Bool)::Bool = nlcache.new_W = val
get_new_W!(nlsolver::AbstractNLSolver)::Bool = get_new_W!(nlsolver.cache)
get_new_W!(nlcache::Union{NLNewtonCache, NLNewtonConstantCache})::Bool = nlcache.new_W
get_new_W!(::AbstractNLSolverCache)::Bool = true

get_W(nlsolver::AbstractNLSolver) = get_W(nlsolver.cache)
get_W(nlcache::Union{NLNewtonCache, NLNewtonConstantCache}) = nlcache.W

set_W_γdt!(nlsolver::AbstractNLSolver, W_γdt) = set_W_γdt!(nlsolver.cache, W_γdt)
function set_W_γdt!(nlcache::Union{NLNewtonCache, NLNewtonConstantCache}, W_γdt)
    nlcache.W_γdt = W_γdt
    W_γdt
end

du_cache(nlsolver::AbstractNLSolver) = du_cache(nlsolver.cache)
du_cache(::AbstractNLSolverCache) = nothing
du_cache(nlcache::Union{NLFunctionalCache, NLAndersonCache, NLNewtonCache}) = (nlcache.k,)

function du_alias_or_new(nlsolver::AbstractNLSolver, rate_prototype)
    _du_cache = du_cache(nlsolver)
    if _du_cache === nothing
        zero(rate_prototype)
    else
        first(_du_cache)
    end
end

mutable struct DAEResidualJacobianWrapper{isAD, F, pType, duType, uType, alphaType,
    gammaType,
    tmpType, uprevType, tType} <: Function
    f::F
    p::pType
    tmp_du::duType
    tmp_u::uType
    α::alphaType
    invγdt::gammaType
    tmp::tmpType
    uprev::uprevType
    t::tType
    function DAEResidualJacobianWrapper(alg, f, p, α, invγdt, tmp, uprev, t)
        isautodiff = alg_autodiff(alg) isa AutoForwardDiff
        if isautodiff
            tmp_du = PreallocationTools.dualcache(uprev)
            tmp_u = PreallocationTools.dualcache(uprev)
        else
            tmp_du = similar(uprev)
            tmp_u = similar(uprev)
        end
        new{isautodiff, typeof(f), typeof(p), typeof(tmp_du), typeof(tmp_u), typeof(α),
            typeof(invγdt), typeof(tmp), typeof(uprev), typeof(t)}(f, p, tmp_du, tmp_u, α,
            invγdt, tmp, uprev, t)
    end
end

is_autodiff(m::DAEResidualJacobianWrapper{isAD}) where {isAD} = isAD

function (m::DAEResidualJacobianWrapper)(out, x)
    if is_autodiff(m)
        tmp_du = PreallocationTools.get_tmp(m.tmp_du, x)
        tmp_u = PreallocationTools.get_tmp(m.tmp_u, x)
    else
        tmp_du = m.tmp_du
        tmp_u = m.tmp_u
    end
    @. tmp_du = (m.α * x + m.tmp) * m.invγdt
    @. tmp_u = x + m.uprev
    m.f(out, tmp_du, tmp_u, m.p, m.t)
end

mutable struct DAEResidualDerivativeWrapper{F, pType, alphaType, gammaType, tmpType,
    uprevType, tType}
    f::F
    p::pType
    α::alphaType
    invγdt::gammaType
    tmp::tmpType
    uprev::uprevType
    t::tType
end

function (m::DAEResidualDerivativeWrapper)(x)
    tmp_du = (m.α * x + m.tmp) * m.invγdt
    tmp_u = x + m.uprev
    m.f(tmp_du, tmp_u, m.p, m.t)
end

DiffEqBase.has_jac(f::DAEResidualJacobianWrapper) = DiffEqBase.has_jac(f.f)
DiffEqBase.has_Wfact(f::DAEResidualJacobianWrapper) = DiffEqBase.has_Wfact(f.f)
DiffEqBase.has_Wfact_t(f::DAEResidualJacobianWrapper) = DiffEqBase.has_Wfact_t(f.f)

DiffEqBase.has_jac(f::DAEResidualDerivativeWrapper) = DiffEqBase.has_jac(f.f)
DiffEqBase.has_Wfact(f::DAEResidualDerivativeWrapper) = DiffEqBase.has_Wfact(f.f)
DiffEqBase.has_Wfact_t(f::DAEResidualDerivativeWrapper) = DiffEqBase.has_Wfact_t(f.f)

function build_nlsolver(alg, u, uprev, p, t, dt, f::F, rate_prototype,
        ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, γ, c,
        iip) where {F, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits,
        tTypeNoUnits, γ, c, 1, iip)
end

function build_nlsolver(alg, u, uprev, p, t, dt, f::F, rate_prototype,
        ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, γ, c, α,
        iip) where {F, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    build_nlsolver(alg, alg.nlsolve, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, α, iip)
end

function build_nlsolver(
        alg, nlalg::Union{NLFunctional, NLAnderson, NLNewton, NonlinearSolveAlg},
        u, uprev, p, t, dt,
        f::F, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        γ, c, α,
        ::Val{true}) where {F, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits}
    #TODO
    #nlalg = DiffEqBase.handle_defaults(alg, nlalg)
    # define unitless type
    uTolType = real(uBottomEltypeNoUnits)
    isdae = alg isa DAEAlgorithm

    # define fields of non-linear solver
    z = zero(u)
    tmp = zero(u)
    ztmp = zero(u)

    # build cache of non-linear solver
    ustep = zero(u)
    tstep = zero(t)
    k = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    dz = zero(u)

    if nlalg isa Union{NLNewton, NonlinearSolveAlg}
        nf = nlsolve_f(f, alg)
        J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))

        # TODO: check if the solver is iterative
        weight = zero(u)

        if islinear(f) || SciMLBase.has_jac(f)
            du1 = rate_prototype
            uf = nothing
            jac_config = nothing
        else
            du1 = zero(rate_prototype)
            if isdae
                uf = DAEResidualJacobianWrapper(alg, f, p, α, inv(γ * dt), tmp, uprev, t)
            else
                uf = build_uf(alg, nf, t, p, Val(true))
            end
            jac_config = build_jac_config(alg, nf, uf, du1, uprev, u, ztmp, dz)
        end
        linprob = LinearProblem(W, _vec(k); u0 = _vec(dz))
        Pl, Pr = wrapprecs(
            alg.precs(W, nothing, u, p, t, nothing, nothing, nothing,
                nothing)...,
            weight, dz)
        linsolve = init(linprob, alg.linsolve,
            alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
            Pl = Pl, Pr = Pr,
            assumptions = LinearSolve.OperatorAssumptions(true))

        tType = typeof(t)
        invγdt = inv(oneunit(t) * one(uTolType))

        if nlalg isa NonlinearSolveAlg
            α = tTypeNoUnits(α)
            dt = tTypeNoUnits(dt)
            if isdae
                nlf = (ztmp, z, p) -> begin
                    tmp, ustep, γ, α, tstep, k, invγdt, _p, dt, f = p
                    _compute_rhs!(tmp, ztmp, ustep, γ, α, tstep, k, invγdt, _p, dt, f, z)[1]
                end
                nlp_params = (tmp, ustep, γ, α, tstep, k, invγdt, p, dt, f)
            else
                nlf = (ztmp, z, p) -> begin
                    tmp, ustep, γ, α, tstep, k, invγdt, method, _p, dt, f = p
                    _compute_rhs!(
                        tmp, ztmp, ustep, γ, α, tstep, k, invγdt, method, _p, dt, f, z)[1]
                end
                nlp_params = (tmp, ustep, γ, α, tstep, k, invγdt, DIRK, p, dt, f)
            end
            prob = NonlinearProblem(NonlinearFunction(nlf), ztmp, nlp_params)
            cache = init(prob, nlalg.alg)
            nlcache = NonlinearSolveCache(ustep, tstep, k, atmp, invγdt, prob, cache)
        else
            nlcache = NLNewtonCache(ustep, tstep, k, atmp, dz, J, W, true,
                true, true, tType(dt), du1, uf, jac_config,
                linsolve, weight, invγdt, tType(nlalg.new_W_dt_cutoff), t)
        end
    elseif nlalg isa NLFunctional
        nlcache = NLFunctionalCache(ustep, tstep, k, atmp, dz)
    elseif nlalg isa NLAnderson
        max_history = min(nlalg.max_history, nlalg.max_iter, length(z))
        Δz₊s = [zero(z) for i in 1:max_history]
        Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
        R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
        γs = Vector{uEltypeNoUnits}(undef, max_history)

        dzold = zero(z)
        z₊old = zero(z)

        nlcache = NLAndersonCache(ustep, tstep, atmp, k, dz, dzold, z₊old, Δz₊s, Q, R, γs,
            0,
            nlalg.aa_start, nlalg.droptol)
    end

    # build non-linear solver
    ηold = one(t)

    NLSolver{true, tTypeNoUnits}(z, tmp, ztmp, γ, c, α, nlalg, nlalg.κ,
        nlalg.fast_convergence_cutoff, ηold, 0, nlalg.max_iter,
        Divergence, nlcache)
end

function build_nlsolver(
        alg, nlalg::Union{NLFunctional, NLAnderson, NLNewton, NonlinearSolveAlg},
        u, uprev, p,
        t, dt,
        f::F, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        γ, c, α,
        ::Val{false}) where {F, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits}
    #TODO
    #nlalg = DiffEqBase.handle_defaults(alg, nlalg)
    # define unitless type
    uTolType = real(uBottomEltypeNoUnits)
    isdae = alg isa DAEAlgorithm

    # define fields of non-linear solver
    z = u
    tmp = u
    ztmp = u

    # build cache of non-linear solver
    tstep = zero(t)

    if nlalg isa Union{NLNewton, NonlinearSolveAlg}
        nf = nlsolve_f(f, alg)
        if isdae
            uf = DAEResidualDerivativeWrapper(f, p, α, inv(γ * dt), tmp, uprev, t)
        else
            uf = build_uf(alg, nf, t, p, Val(false))
        end

        tType = typeof(t)
        invγdt = inv(oneunit(t) * one(uTolType))

        J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
        if nlalg isa NonlinearSolveAlg
            α = tTypeNoUnits(α)
            dt = tTypeNoUnits(dt)
            if isdae
                nlf = (z, p) -> begin
                    tmp, α, tstep, invγdt, _p, dt, uprev, f = p
                    _compute_rhs(tmp, α, tstep, invγdt, p, dt, uprev, f, z)[1]
                end
                nlp_params = (tmp, α, tstep, invγdt, _p, dt, uprev, f)
            else
                nlf = (z, p) -> begin
                    tmp, γ, α, tstep, invγdt, method, _p, dt, f = p
                    _compute_rhs(tmp, γ, α, tstep, invγdt, method, _p, dt, f, z)[1]
                end
                nlp_params = (tmp, γ, α, tstep, invγdt, DIRK, p, dt, f)
            end
            prob = NonlinearProblem(NonlinearFunction(nlf), copy(ztmp), nlp_params)
            cache = init(prob, nlalg.alg)
            nlcache = NonlinearSolveCache(
                nothing, tstep, nothing, nothing, invγdt, prob, cache)
        else
            nlcache = NLNewtonConstantCache(tstep, J, W, true, true, true, tType(dt), uf,
                invγdt, tType(nlalg.new_W_dt_cutoff), t)
        end
    elseif nlalg isa NLFunctional
        nlcache = NLFunctionalConstantCache(tstep)
    elseif nlalg isa NLAnderson
        max_history = min(nlalg.max_history, nlalg.max_iter, length(z))
        Δz₊s = Vector{typeof(z)}(undef, max_history)
        Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
        R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
        γs = Vector{uEltypeNoUnits}(undef, max_history)

        dz = u
        dzold = u
        z₊old = u

        nlcache = NLAndersonConstantCache(tstep, dz, dzold, z₊old, Δz₊s, Q, R, γs, 0,
            nlalg.aa_start, nlalg.droptol)
    end

    # build non-linear solver
    ηold = one(tTypeNoUnits)
    NLSolver{false, tTypeNoUnits}(z, tmp, ztmp, γ, c, α, nlalg, nlalg.κ,
        nlalg.fast_convergence_cutoff, ηold, 0, nlalg.max_iter,
        Divergence,
        nlcache)
end

## Anderson acceleration

"""
    anderson(z, cache)

Return the next iterate of the fixed-point iteration `z = g(z)` by performing Anderson
acceleration based on the current iterate `z` and the settings and history in the `cache`.
"""
@muladd function anderson(z, cache)
    @unpack dz, Δz₊s, z₊old, dzold, R, Q, γs, history, droptol = cache

    # increase size of history
    history += 1

    # remove oldest history if maximum size is exceeded
    max_history = length(Δz₊s)
    if history > max_history
        # circularly shift differences of z₊
        for i in 1:(max_history - 1)
            Δz₊s[i] = Δz₊s[i + 1]
        end

        # delete left-most column of QR decomposition
        qrdelete!(Q, R, max_history)

        # update size of history
        history = max_history
    end

    # update history of differences of z₊
    Δz₊s[history] = @.. broadcast=false z-z₊old

    # replace/add difference of residuals as right-most column to QR decomposition
    qradd!(Q, R, _vec(dz .- dzold), history)

    # update cached values
    cache.dzold = dz
    cache.z₊old = z

    # define current Q and R matrices
    Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

    # check condition (TODO: incremental estimation)
    if droptol !== nothing
        while cond(R) > droptol && history > 1
            qrdelete!(Q, R, history)
            history -= 1
            Qcur, Rcur = view(Q, :, 1:history),
            UpperTriangular(view(R, 1:history, 1:history))
        end
    end

    # save updated history
    cache.history = history

    # solve least squares problem
    γscur = view(γs, 1:history)
    ldiv!(Rcur, mul!(γscur, Qcur', _vec(dz)))

    # update next iterate
    for i in 1:history
        z = @.. broadcast=false z-γs[i] * Δz₊s[i]
    end

    z
end

"""
    anderson!(z, cache)

Update the current iterate `z` of the fixed-point iteration `z = g(z)` in-place
by performing Anderson acceleration based on the settings and history in the `cache`.
"""
@muladd function anderson!(z, cache)
    @unpack dz, z₊old, dzold, Δz₊s, γs, R, Q, history, droptol = cache

    # increase size of history
    history += 1

    # remove oldest history if maximum size is exceeded
    max_history = length(Δz₊s)
    if history > max_history
        # circularly shift differences of z₊
        ptr = Δz₊s[1]
        for i in 1:(max_history - 1)
            Δz₊s[i] = Δz₊s[i + 1]
        end
        Δz₊s[max_history] = ptr

        # delete left-most column of QR decomposition
        qrdelete!(Q, R, max_history)

        # update size of history
        history = max_history
    end

    # update history of differences of z₊
    @.. broadcast=false Δz₊s[history]=z - z₊old

    # replace/add difference of residuals as right-most column to QR decomposition
    @.. broadcast=false dzold=dz - dzold
    qradd!(Q, R, _vec(dzold), history)

    # update cached values
    @.. broadcast=false dzold=dz
    @.. broadcast=false z₊old=z

    # define current Q and R matrices
    Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

    # check condition (TODO: incremental estimation)
    if droptol !== nothing
        while cond(R) > droptol && history > 1
            qrdelete!(Q, R, history)
            history -= 1
            Qcur, Rcur = view(Q, :, 1:history),
            UpperTriangular(view(R, 1:history, 1:history))
        end
    end

    # save updated history
    cache.history = history

    # solve least squares problem
    γscur = view(γs, 1:history)
    ldiv!(Rcur, mul!(γscur, Qcur', _vec(dz)))

    # update next iterate
    for i in 1:history
        @.. broadcast=false z=z - γs[i] * Δz₊s[i]
    end

    nothing
end

## resize

function resize_nlsolver!(integrator::DiffEqBase.DEIntegrator, i::Int)
    isdefined(integrator.cache, :nlsolver) || return

    @unpack nlsolver = integrator.cache

    if nlsolver isa AbstractArray
        for idx in eachindex(nlsolver)
            resize!(nlsolver[idx], integrator, i)
        end
    else
        resize!(nlsolver, integrator, i)
    end

    nlsolver.alg isa NLNewton && resize!(nlsolver.cache.linsolve, i)

    # make it reset everything since the caches changed size!
    nlsolver.cache.firstcall = true

    nothing
end

function Base.resize!(nlsolver::AbstractNLSolver, integrator, i::Int)
    resize!(nlsolver.z, i)
    resize!(nlsolver.tmp, i)
    resize!(nlsolver.ztmp, i)

    resize!(nlsolver.cache, nlsolver, integrator, i)
end

## default: dispatch only on the cache
function Base.resize!(cache::AbstractNLSolverCache, nlsolver, integrator, i::Int)
    Base.resize!(cache, i)
end

"""
    qrdelete!(Q, R, k)

Delete the left-most column of F = Q[:, 1:k] * R[1:k, 1:k] by updating Q and R.
Only Q[:, 1:(k-1)] and R[1:(k-1), 1:(k-1)] are valid on exit.
"""
function qrdelete!(Q::AbstractMatrix, R::AbstractMatrix, k::Int)
    n, m = size(Q)
    m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
    1 ≤ k ≤ m || throw(ArgumentError())

    # apply Givens rotations
    for i in 2:k
        g = first(givens(R, i - 1, i, i))
        lmul!(g, R)
        rmul!(Q, g')
    end

    # move columns of R
    @inbounds for j in 1:(k - 1)
        for i in 1:(k - 1)
            R[i, j] = R[i, j + 1]
        end
    end

    Q, R
end

"""
    qradd!(Q, R, v, k)

Replace the right-most column of F = Q[:, 1:k] * R[1:k, 1:k] with v by updating Q and R.
This implementation modifies vector v as well. Only Q[:, 1:k] and R[1:k, 1:k] are valid on
exit.
"""
function qradd!(Q::AbstractMatrix, R::AbstractMatrix, v::AbstractVector, k::Int)
    n, m = size(Q)
    n == length(v) || throw(DimensionMismatch())
    m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
    1 ≤ k ≤ m || throw(ArgumentError())

    @inbounds for i in 1:(k - 1)
        q = view(Q, :, i)
        r = dot(q, v)

        R[i, k] = r
        axpy!(-r, q, v)
    end

    @inbounds begin
        d = norm(v)
        R[k, k] = d
        @.. broadcast=false @view(Q[:, k])=v / d
    end

    Q, R
end

function qradd!(Q::AbstractMatrix, R::AbstractMatrix, v::Number, k::Int)
    1 == LinearAlgebra.checksquare(Q) == LinearAlgebra.checksquare(R) ||
        throw(DimensionMismatch())
    k == 1 || throw(ArgumentError())

    R[1, 1] = abs(v)
    Q[1, 1] = one(v)

    Q, R
end
