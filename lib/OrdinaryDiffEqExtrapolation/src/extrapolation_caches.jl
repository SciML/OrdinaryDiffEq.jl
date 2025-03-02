abstract type ExtrapolationMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::ExtrapolationMutableCache, u) = (cache.fsalfirst, cache.k)

@cache mutable struct AitkenNevilleCache{
    uType,
    rateType,
    arrayType,
    dtType,
    uNoUnitsType
} <: ExtrapolationMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
    utilde::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    dtpropose::dtType
    T::arrayType
    cur_order::Int
    work::dtType
    A::Int
    step_no::Int
    u_tmps::Array{uType, 1}
    k_tmps::Array{rateType, 1}
end

@cache mutable struct AitkenNevilleConstantCache{dtType, arrayType} <:
                      OrdinaryDiffEqConstantCache
    dtpropose::dtType
    T::arrayType
    cur_order::Int
    work::dtType
    A::Int
    step_no::Int
end

function alg_cache(alg::AitkenNeville, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    utilde = zero(u)
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    cur_order = max(alg.init_order, alg.min_order)
    dtpropose = zero(dt)
    T = Array{typeof(u), 2}(undef, alg.max_order, alg.max_order)
    # Array of arrays of length equal to number of threads to store intermediate
    # values of u and k. [Thread Safety]
    u_tmps = Array{typeof(u), 1}(undef, Threads.nthreads())
    k_tmps = Array{typeof(k), 1}(undef, Threads.nthreads())
    # Initialize each element of u_tmps and k_tmps to different instance of
    # zeros array similar to u and k respectively
    for i in 1:Threads.nthreads()
        u_tmps[i] = zero(u)
        k_tmps[i] = zero(rate_prototype)
    end
    # Initialize lower triangle of T to different instance of zeros array similar to u
    for i in 1:(alg.max_order)
        for j in 1:i
            T[i, j] = zero(u)
        end
    end
    work = zero(dt)
    A = one(Int)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    step_no = zero(Int)
    AitkenNevilleCache(u, uprev, tmp, k, utilde, atmp, fsalfirst, dtpropose, T, cur_order,
        work, A, step_no, u_tmps, k_tmps)
end

function alg_cache(alg::AitkenNeville, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dtpropose = zero(dt)
    cur_order = max(alg.init_order, alg.min_order)
    T = Array{typeof(u), 2}(undef, alg.max_order, alg.max_order)
    @.. broadcast=false T=u
    work = zero(dt)
    A = one(Int)
    step_no = zero(Int)
    AitkenNevilleConstantCache(dtpropose, T, cur_order, work, A, step_no)
end

@cache mutable struct ImplicitEulerExtrapolationCache{uType, rateType, QType, arrayType,
    dtType, JType, WType, F, JCType,
    GCType, uNoUnitsType, TFType, UFType,
    sequenceType} <:
                      ExtrapolationMutableCache
    uprev::uType
    u_tmps::Array{uType, 1}
    u_tmps2::Array{uType, 1}
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    k_tmps::Array{rateType, 1}
    dtpropose::dtType
    T::arrayType
    A::Int
    step_no::Int
    du1::rateType
    du2::rateType
    J::JType
    W::WType
    tf::TFType
    uf::UFType
    linsolve_tmps::Array{rateType, 1}
    linsolve::Array{F, 1}
    jac_config::JCType
    grad_config::GCType
    sequence::sequenceType #support for different sequences
    stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n - 1)

    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
    n_curr::Int # Storage for the current extrapolation order
    n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter
    sigma::Rational{Int} # Parameter for order selection
    res::uNoUnitsType # Storage for the scaled residual of u and utilde

    #Stepsizing caches
    work::Array{QType, 1}
    dt_new::Array{QType, 1}

    # Values to check overflow in T1 computation
    diff1::Array{uType, 1}
    diff2::Array{uType, 1}
end
get_fsalfirstlast(cache::ImplicitEulerExtrapolationCache, u) = (zero(u), zero(u))

@cache mutable struct ImplicitEulerExtrapolationConstantCache{QType, dtType, arrayType, TF,
    UF, sequenceType} <:
                      OrdinaryDiffEqConstantCache
    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n)
    dtpropose::dtType
    T::arrayType
    n_curr::Int
    n_old::Int
    A::Int
    step_no::Int
    sigma::Rational{Int}

    tf::TF
    uf::UF

    sequence::sequenceType #support for different sequences
    stage_number::Vector{Int}

    #Stepsizing caches
    work::Array{QType, 1}
    dt_new::Array{QType, 1}
end

function alg_cache(alg::ImplicitEulerExtrapolation, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dtpropose = zero(dt)
    #cur_order = max(alg.init_order, alg.min_order)
    QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl
    Q = fill(zero(QType), alg.max_order + 1)
    n_curr = alg.init_order
    n_old = alg.init_order
    T = Array{typeof(u), 2}(undef, alg.max_order + 1, alg.max_order + 1)
    for i in 1:(alg.max_order + 1)
        for j in 1:i
            T[i, j] = zero(u)
        end
    end
    A = one(Int)
    step_no = zero(Int)
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    sequence = generate_sequence(constvalue(uBottomEltypeNoUnits), alg)
    stage_number = Vector{Int}(undef, alg.max_order + 1)

    for n in 1:length(stage_number)
        s = zero(eltype(sequence))
        for i in 1:n
            s += sequence[i]
        end
        stage_number[n] = 2 * Int(s) - n + 7
    end
    sigma = 9 // 10
    work = fill(zero(eltype(Q)), alg.max_order + 1)
    dt_new = fill(zero(eltype(Q)), alg.max_order + 1)
    ImplicitEulerExtrapolationConstantCache(Q, dtpropose, T, n_curr, n_old, A, step_no,
        sigma, tf, uf, sequence, stage_number, work,
        dt_new)
end

function alg_cache(alg::ImplicitEulerExtrapolation, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    u_tmp = zero(u)
    u_tmps = Array{typeof(u_tmp), 1}(undef, Threads.nthreads())

    u_tmps[1] = u_tmp
    for i in 2:Threads.nthreads()
        u_tmps[i] = zero(u_tmp)
    end

    u_tmps2 = Array{typeof(u_tmp), 1}(undef, Threads.nthreads())

    for i in 1:Threads.nthreads()
        u_tmps2[i] = zero(u_tmp)
    end

    utilde = zero(u)
    tmp = zero(u)
    k_tmp = zero(rate_prototype)
    k_tmps = Array{typeof(k_tmp), 1}(undef, Threads.nthreads())

    k_tmps[1] = k_tmp
    for i in 2:Threads.nthreads()
        k_tmps[i] = zero(rate_prototype)
    end

    #cur_order = max(alg.init_order, alg.min_order)
    dtpropose = zero(dt)
    T = Array{typeof(u), 2}(undef, alg.max_order + 1, alg.max_order + 1)
    # Initialize lower triangle of T to different instance of zeros array similar to u
    for i in 1:(alg.max_order + 1)
        for j in 1:i
            T[i, j] = zero(u)
        end
    end
    A = one(Int)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    step_no = zero(Int)

    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)

    if DiffEqBase.has_jac(f) && !DiffEqBase.has_Wfact(f) && f.jac_prototype !== nothing
        W_el = WOperator(f, dt, true)
        J = nothing # is J = W.J better?
    else
        J = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
        W_el = zero(J)
    end

    W = Array{typeof(W_el), 1}(undef, Threads.nthreads())
    W[1] = W_el
    for i in 2:Threads.nthreads()
        if W_el isa WOperator
            W[i] = WOperator(f, dt, true)
        else
            W[i] = zero(W_el)
        end
    end

    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)
    linsolve_tmp = zero(rate_prototype)
    linsolve_tmps = Array{typeof(linsolve_tmp), 1}(undef, Threads.nthreads())

    for i in 1:Threads.nthreads()
        linsolve_tmps[i] = zero(rate_prototype)
    end

    linprob = LinearProblem(W[1], _vec(linsolve_tmps[1]); u0 = _vec(k_tmps[1]))
    linsolve1 = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true))
    #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
    #Pr = Diagonal(_vec(weight)))

    linsolve = Array{typeof(linsolve1), 1}(undef, Threads.nthreads())
    linsolve[1] = linsolve1
    for i in 2:Threads.nthreads()
        linprob = LinearProblem(W[i], _vec(linsolve_tmps[i]); u0 = _vec(k_tmps[i]))
        linsolve[i] = init(linprob, alg.linsolve,
            alias = LinearAliasSpecifier(alias_A = true, alias_b = true))
        #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
        #Pr = Diagonal(_vec(weight)))
    end

    res = uEltypeNoUnits.(zero(u))
    grad_config = build_grad_config(alg, f, tf, du1, t)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, du1, du2)
    sequence = generate_sequence(constvalue(uBottomEltypeNoUnits), alg)
    cc = alg_cache(alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(false))
    diff1 = Array{typeof(u), 1}(undef, Threads.nthreads())
    diff2 = Array{typeof(u), 1}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        diff1[i] = zero(u)
        diff2[i] = zero(u)
    end
    ImplicitEulerExtrapolationCache(uprev, u_tmps, u_tmps2, utilde, tmp, atmp, k_tmps,
        dtpropose, T, A, step_no,
        du1, du2, J, W, tf, uf, linsolve_tmps, linsolve,
        jac_config, grad_config, sequence, cc.stage_number,
        cc.Q, cc.n_curr, cc.n_old, cc.sigma, res, cc.work,
        cc.dt_new, diff1, diff2)
end

struct extrapolation_coefficients{T1, T2, T3}
    # This structure is used by the caches of the algorithms
    # ExtrapolationMidpointDeuflhard() and  ExtrapolationMidpointHairerWanner().
    # It contains the constant coefficients used to extrapolate the internal discretisations
    # in their perfom_step! function and some additional constant data.

    subdividing_sequence::T1  # subdividing_sequence[n] is used for the (n -1)th internal discretisation

    # Weights and Scaling factors for extrapolation operators
    extrapolation_weights::T2
    extrapolation_scalars::T3

    # Weights and scaling factors for internal extrapolation operators (used for error estimate)
    extrapolation_weights_2::T2
    extrapolation_scalars_2::T3
end

function create_extrapolation_coefficients(T,
        alg::Union{ExtrapolationMidpointDeuflhard,
            ExtrapolationMidpointHairerWanner,
            ImplicitDeuflhardExtrapolation,
            ImplicitHairerWannerExtrapolation})
    # Compute and return extrapolation_coefficients

    @unpack min_order, init_order, max_order, sequence = alg

    # Initialize subdividing_sequence:
    if sequence == :harmonic
        subdividing_sequence = BigInt.(1:(max_order + 1))
    elseif sequence == :romberg
        subdividing_sequence = BigInt(2) .^ (0:max_order)
    else # sequence == :bulirsch
        subdividing_sequence = [n == 0 ? BigInt(1) :
                                (isodd(n) ? BigInt(2)^((n + 1) ÷ 2) :
                                 3 * BigInt(2)^(n ÷ 2 - 1)) for n in 0:max_order]
    end

    # Compute nodes corresponding to subdividing_sequence
    nodes = BigInt(1) .// subdividing_sequence .^ 2

    # Compute barycentric weights for internal extrapolation operators
    extrapolation_weights_2 = zeros(Rational{BigInt}, max_order, max_order)
    extrapolation_weights_2[1, :] = ones(Rational{BigInt}, 1, max_order)
    for n in 2:max_order
        distance = nodes[2:n] .- nodes[n + 1]
        extrapolation_weights_2[1:(n - 1), n] = extrapolation_weights_2[1:(n - 1),
            n - 1] .// distance
        extrapolation_weights_2[n, n] = 1 // prod(-distance)
    end

    # Compute barycentric weights for extrapolation operators
    extrapolation_weights = zeros(Rational{BigInt}, max_order + 1, max_order + 1)
    for n in 1:max_order
        extrapolation_weights[n + 1, (n + 1):(max_order + 1)] = extrapolation_weights_2[n,
            n:max_order] //
                                                                (nodes[n + 1] - nodes[1])
        extrapolation_weights[1, n] = 1 // prod(nodes[1] .- nodes[2:n])
    end
    extrapolation_weights[1, max_order + 1] = 1 //
                                              prod(nodes[1] .- nodes[2:(max_order + 1)])

    # Rescale barycentric weights to obtain weights of 1. Barycentric Formula
    for m in 1:(max_order + 1)
        extrapolation_weights[1:m, m] = -extrapolation_weights[1:m, m] .// nodes[1:m]
        if 2 <= m
            extrapolation_weights_2[1:(m - 1), m - 1] = -extrapolation_weights_2[1:(m - 1),
                m - 1] .//
                                                        nodes[2:m]
        end
    end

    # Compute scaling factors for internal extrapolation operators
    extrapolation_scalars_2 = ones(Rational{BigInt}, max_order)
    extrapolation_scalars_2[1] = -nodes[2]
    for n in 1:(max_order - 1)
        extrapolation_scalars_2[n + 1] = -extrapolation_scalars_2[n] * nodes[n + 2]
    end

    # Compute scaling factors for extrapolation operators
    extrapolation_scalars = -nodes[1] * [BigInt(1); extrapolation_scalars_2]

    # Initialize and return extrapolation_coefficients
    extrapolation_coefficients(Int.(subdividing_sequence),
        T.(extrapolation_weights), T.(extrapolation_scalars),
        T.(extrapolation_weights_2), T.(extrapolation_scalars_2))
end

function create_extrapolation_coefficients(T, alg::ImplicitEulerBarycentricExtrapolation)
    # Compute and return extrapolation_coefficients

    @unpack min_order, init_order, max_order, sequence = alg

    # Initialize subdividing_sequence:
    if sequence == :harmonic
        subdividing_sequence = BigInt.(1:(max_order + 1))
    elseif sequence == :romberg
        subdividing_sequence = BigInt(2) .^ (0:max_order)
    else # sequence == :bulirsch
        subdividing_sequence = [n == 0 ? BigInt(1) :
                                (isodd(n) ? BigInt(2)^((n + 1) ÷ 2) :
                                 3 * BigInt(2)^(n ÷ 2 - 1)) for n in 0:max_order]
    end

    # Compute nodes corresponding to subdividing_sequence
    nodes = BigInt(1) .// subdividing_sequence

    # Compute barycentric weights for internal extrapolation operators
    extrapolation_weights_2 = zeros(Rational{BigInt}, max_order, max_order)
    extrapolation_weights_2[1, :] = ones(Rational{BigInt}, 1, max_order)
    for n in 2:max_order
        distance = nodes[2:n] .- nodes[n + 1]
        extrapolation_weights_2[1:(n - 1), n] = extrapolation_weights_2[1:(n - 1),
            n - 1] .// distance
        extrapolation_weights_2[n, n] = 1 // prod(-distance)
    end

    # Compute barycentric weights for extrapolation operators
    extrapolation_weights = zeros(Rational{BigInt}, max_order + 1, max_order + 1)
    for n in 1:max_order
        extrapolation_weights[n + 1, (n + 1):(max_order + 1)] = extrapolation_weights_2[n,
            n:max_order] //
                                                                (nodes[n + 1] - nodes[1])
        extrapolation_weights[1, n] = 1 // prod(nodes[1] .- nodes[2:n])
    end
    extrapolation_weights[1, max_order + 1] = 1 //
                                              prod(nodes[1] .- nodes[2:(max_order + 1)])

    # Rescale barycentric weights to obtain weights of 1. Barycentric Formula
    for m in 1:(max_order + 1)
        extrapolation_weights[1:m, m] = -extrapolation_weights[1:m, m] .// nodes[1:m]
        if 2 <= m
            extrapolation_weights_2[1:(m - 1), m - 1] = -extrapolation_weights_2[1:(m - 1),
                m - 1] .//
                                                        nodes[2:m]
        end
    end

    # Compute scaling factors for internal extrapolation operators
    extrapolation_scalars_2 = ones(Rational{BigInt}, max_order)
    extrapolation_scalars_2[1] = -nodes[2]
    for n in 1:(max_order - 1)
        extrapolation_scalars_2[n + 1] = -extrapolation_scalars_2[n] * nodes[n + 2]
    end

    # Compute scaling factors for extrapolation operators
    extrapolation_scalars = -nodes[1] * [BigInt(1); extrapolation_scalars_2]

    # Initialize and return extrapolation_coefficients
    extrapolation_coefficients(Int.(subdividing_sequence),
        T.(extrapolation_weights), T.(extrapolation_scalars),
        T.(extrapolation_weights_2), T.(extrapolation_scalars_2))
end

function create_extrapolation_coefficients(T::Type{<:CompiledFloats},
        alg::ImplicitEulerBarycentricExtrapolation)
    # Compute and return extrapolation_coefficients

    @unpack min_order, init_order, max_order, sequence = alg

    max_order > 15 &&
        error("max_order > 15 not allowed for Float32 or Float64 with this algorithm. That's a bad idea.")

    # Initialize subdividing_sequence:
    if sequence == :harmonic
        subdividing_sequence = [
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21
        ]
        extrapolation_weights = T[-1.0 -2.0 -3.0 -4.0 -5.0 -6.0 -7.0 -8.0 -9.0 -10.0 -11.0 -12.0 -13.0 -14.0 -15.0 -16.0;
                                  0.0 4.0 24.0 96.0 320.0 960.0 2688.0 7168.0 18432.0 46080.0 112640.0 270336.0 638976.0 1.490944e6 3.44064e6 7.86432e6;
                                  0.0 0.0 -27.0 -324.0 -2430.0 -14580.0 -76545.0 -367416.0 -1.653372e6 -7.08588e6 -2.9229255e7 -1.1691702e8 -4.55976378e8 -1.741000716e9 -6.528752685e9 -2.410616376e10;
                                  0.0 0.0 0.0 256.0 5120.0 61440.0 573440.0 4.58752e6 3.3030144e7 2.2020096e8 1.38412032e9 8.30472192e9 4.798283776e10 2.68703891456e11 1.46565758976e12 7.81684047872e12;
                                  0.0 0.0 0.0 0.0 -3125.0 -93750.0 -1.640625e6 -2.1875e7 -2.4609375e8 -2.4609375e9 -2.255859375e10 -1.93359375e11 -1.571044921875e12 -1.221923828125e13 -9.1644287109375e13 -6.6650390625e14;
                                  0.0 0.0 0.0 0.0 0.0 46656.0 1.959552e6 4.7029248e7 8.46526464e8 1.269789696e10 1.67612239872e11 2.011346878464e12 2.2412150931456e13 2.35327584780288e14 2.35327584780288e15 2.259144813890765e16;
                                  0.0 0.0 0.0 0.0 0.0 0.0 -823543.0 -4.6118408e7 -1.452729852e9 -3.389702988e10 -6.5251782519e11 -1.0962299463192e13 -1.66261541858412e14 -2.327661586017768e15 -3.0550558316483204e16 -3.8018472571623546e17;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.6777216e7 1.207959552e9 4.831838208e10 1.41733920768e12 3.401614098432e13 7.07535732473856e14 1.3207333672845312e16 2.2641143439163392e17 3.6225829502661427e18;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -3.87420489e8 -3.486784401e10 -1.725958278495e12 -6.213449802582e13 -1.817434067255235e15 -4.579933849483192e16 -1.0304851161337183e18 -2.119855096046506e19;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0e10 1.1e12 6.6e13 2.86e15 1.001e17 3.003e18 8.008e19;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.85311670611e11 -3.7661140520652e13 -2.692771547226618e15 -1.3822893942429973e17 -5.701943751252363e18 -2.007084200440832e20;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 8.916100448256e12 1.390911669927936e15 1.1683658027394662e17 7.010194816436797e18 3.364893511889663e20;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -3.02875106592253e14 -5.512326939979005e16 -5.37451876647953e18 -3.726333011425807e20;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.1112006825558016e16 2.3335214333671834e18 2.6135440053712454e20;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -4.378938903808594e17 -1.0509453369140625e20;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.8446744073709552e19]
        extrapolation_weights_2 = T[-2.0 -12.0 -48.0 -160.0 -480.0 -1344.0 -3584.0 -9216.0 -23040.0 -56320.0 -135168.0 -319488.0 -745472.0 -1.72032e6 -3.93216e6;
                                    0.0 18.0 216.0 1620.0 9720.0 51030.0 244944.0 1.102248e6 4.72392e6 1.948617e7 7.794468e7 3.03984252e8 1.160667144e9 4.35250179e9 1.607077584e10;
                                    0.0 0.0 -192.0 -3840.0 -46080.0 -430080.0 -3.44064e6 -2.4772608e7 -1.6515072e8 -1.03809024e9 -6.22854144e9 -3.598712832e10 -2.01527918592e11 -1.09924319232e12 -5.86263035904e12;
                                    0.0 0.0 0.0 2500.0 75000.0 1.3125e6 1.75e7 1.96875e8 1.96875e9 1.8046875e10 1.546875e11 1.2568359375e12 9.775390625e12 7.33154296875e13 5.33203125e14;
                                    0.0 0.0 0.0 0.0 -38880.0 -1.63296e6 -3.919104e7 -7.0543872e8 -1.05815808e10 -1.3967686656e11 -1.67612239872e12 -1.867679244288e13 -1.9610632065024e14 -1.9610632065024e15 -1.882620678242304e16;
                                    0.0 0.0 0.0 0.0 0.0 705894.0 3.9530064e7 1.245197016e9 2.905459704e10 5.5930099302e11 9.396256682736e12 1.42509893021496e14 1.995138502300944e15 2.618619284269989e16 3.2587262204248755e17;
                                    0.0 0.0 0.0 0.0 0.0 0.0 -1.4680064e7 -1.056964608e9 -4.227858432e10 -1.24017180672e12 -2.976412336128e13 -6.19093765914624e14 -1.1556416963739648e16 -1.9811000509267968e17 -3.169760081482875e18;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.44373768e8 3.099363912e10 1.53418513644e12 5.523066491184e13 1.61549694867132e15 4.071052310651726e16 9.159867698966385e17 1.884315640930228e19;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -9.0e9 -9.9e11 -5.94e13 -2.574e15 -9.009e16 -2.7027e18 -7.2072e19;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.5937424601e11 3.423740047332e13 2.44797413384238e15 1.2566267220390883e17 5.183585228411239e18 1.8246220004007562e20;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -8.173092077568e12 -1.275002364100608e15 -1.0710019858445107e17 -6.426011915067064e18 -3.084485719232191e20;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.79577021469772e14 5.0883017907498504e16 4.961094245981104e18 3.439692010546899e20;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0318292052303872e16 -2.166841330983813e18 -2.4268622907018707e20;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.0870096435546874e17 9.80882314453125e19;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.7293822569102705e19]
        extrapolation_scalars = T[-1.0, 0.5, -0.16666666666666666, 0.041666666666666664,
            -0.008333333333333333, 0.001388888888888889,
            -0.0001984126984126984, 2.48015873015873e-5,
            -2.7557319223985893e-6, 2.755731922398589e-7,
            -2.505210838544172e-8, 2.08767569878681e-9,
            -1.6059043836821613e-10, 1.1470745597729725e-11,
            -7.647163731819816e-13, 4.779477332387385e-14]
        extrapolation_scalars_2 = T[-0.5, 0.16666666666666666, -0.041666666666666664,
            0.008333333333333333, -0.001388888888888889,
            0.0001984126984126984, -2.48015873015873e-5,
            2.7557319223985893e-6, -2.755731922398589e-7,
            2.505210838544172e-8, -2.08767569878681e-9,
            1.6059043836821613e-10, -1.1470745597729725e-11,
            7.647163731819816e-13, -4.779477332387385e-14]
    elseif sequence == :romberg
        subdividing_sequence = [
            1,
            2,
            4,
            8,
            16,
            32,
            64,
            128,
            256,
            512,
            1024,
            2048,
            4096,
            8192,
            16384,
            32768
        ]
        extrapolation_weights = T[-1.0 -2.0 -2.6666666666666665 -3.0476190476190474 -3.250793650793651 -3.355657962109575 -3.4089223742065524 -3.4357642826648718 -3.4492378680870868 -3.4559878443455743 -3.4593661315834487 -3.461056100382464 -3.4619012911273677 -3.462323938092467 -3.4625352744739657 -3.4626409458895506;
                                  0.0 4.0 16.0 42.666666666666664 97.52380952380952 208.05079365079365 429.5242191500256 872.6841277968774 1759.1113127244143 3532.019576921177 7077.863105219736 14169.563674965806 28352.971574333144 56719.79075383079 113453.43080341395 226920.71174792582;
                                  0.0 0.0 -21.333333333333332 -170.66666666666666 -910.2222222222222 -4161.015873015873 -17753.667724867726 -73305.4667349377 -297876.1822880008 -1.2008866561532002e6 -4.822384062356381e6 -1.9327284852653358e7 -7.738471041687992e7 -3.096900575159161e8 -1.2390627356143517e9 -4.956856027421692e9;
                                  0.0 0.0 0.0 195.04761904761904 3120.7619047619046 33288.12698412698 304348.589569161 2.597107964323507e6 2.144708512473606e7 1.7430012037880734e8 1.4053804981724308e9 1.1287134353949562e10 9.047378143596361e10 7.244977688400918e11 5.798813602675166e12 4.6401837394984086e13;
                                  0.0 0.0 0.0 0.0 -3328.8126984126984 -106522.00634920635 -2.272469468783069e6 -4.1553727429176114e7 -7.091836147912724e8 -1.1712968089455853e10 -1.9038221148575864e11 -3.070100544274281e12 -4.9314242468029234e13 -7.905719653583081e14 -1.2661516207654468e16 -2.026832257412905e17;
                                  0.0 0.0 0.0 0.0 0.0 109958.20010240655 7.037324806554019e6 3.002591917463048e8 1.0980907583864862e10 3.748149788625873e11 1.2380985108235143e13 4.024802778042154e14 1.2980781243197372e16 4.17013960565776e17 1.3370561115283118e19 4.282761941599191e20;
                                  0.0 0.0 0.0 0.0 0.0 0.0 -7.14902837491202e6 -9.150756319887385e8 -7.808645392970569e10 -5.711466344572759e12 -3.899027691228337e14 -2.5758737779469784e16 -1.6747268245191785e18 -1.0802647359418198e20 -6.940806836733636e21 -4.45080936254575e23;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 9.222809519256577e8 2.3610392369296838e11 4.029506964359994e13 5.894593045006619e15 8.04808437078237e17 1.0633855994427287e20 1.3827388620055293e22 1.7838420090628812e24 2.2922719589400975e26;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.370298214329408e11 -1.2135926857366569e14 -4.142396367314456e16 -1.2119468228942864e19 -3.309422791049998e21 -8.745416614284383e23 -2.274363584260878e26 -5.86821621488665e28;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.215967622499351e14 1.2451508454393354e17 8.500229771532529e19 4.97384873488532e22 2.7163845890787027e25 1.435653067982757e28 7.467219005025235e30;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.2463680016909867e17 -2.5525616674631408e20 -3.4850975299763416e23 -4.0785598522237415e26 -4.4548749745889186e29 -4.708946553784829e32;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.5538086443402602e20 1.0460400207217706e24 2.856386616584248e27 6.685576903730903e30 1.4604865598763616e34;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.046295463950274e24 -8.571252440680644e27 -4.68104666627039e31 -2.1912648165764018e35;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 8.572298863881802e27 1.4044854458583944e32 1.5340726363295958e36;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.4045711740794687e32 -4.602498823223603e36;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.602639284627553e36]
        extrapolation_weights_2 = T[-2.0 -8.0 -21.333333333333332 -48.76190476190476 -104.02539682539683 -214.7621095750128 -436.3420638984387 -879.5556563622072 -1766.0097884605884 -3538.931552609868 -7084.781837482903 -14176.485787166572 -28359.895376915396 -56726.71540170698 -113460.35587396291;
                                    0.0 16.0 128.0 682.6666666666666 3120.7619047619046 13315.250793650794 54979.10005120328 223407.13671600062 900664.9921149001 3.616788046767285e6 1.449546363949002e7 5.803853281265994e7 2.322675431369371e8 9.292970517107637e8 3.7176420205662684e9;
                                    0.0 0.0 -170.66666666666666 -2730.6666666666665 -29127.11111111111 -266305.01587301586 -2.272469468783069e6 -1.876619948414405e7 -1.5251260533145642e8 -1.229707935900877e9 -9.876242559705868e9 -7.916455875646815e10 -6.339355477350803e11 -5.07396190234077e12 -4.060160772061108e13;
                                    0.0 0.0 0.0 3120.7619047619046 99864.38095238095 2.130440126984127e6 3.895661946485261e7 6.648596388668178e8 1.0980907583864862e10 1.784833232678987e11 2.878219260257138e12 4.623210231377741e13 7.411612175234139e14 1.1870171444676064e16 1.9001552413245984e17;
                                    0.0 0.0 0.0 0.0 -106522.00634920635 -6.817408406349206e6 -2.908760920042328e8 -1.0637754221869085e10 -3.631020107731315e11 -1.1994079323602793e13 -3.899027691228337e14 -1.2575131829347454e16 -4.039822742980955e17 -1.295273108043052e19 -4.148925630924216e20;
                                    0.0 0.0 0.0 0.0 0.0 7.037324806554019e6 9.007775752389145e8 7.686635308705403e10 5.62222468293881e12 3.838105383552894e14 2.5356257501665572e16 1.6485592178860662e18 1.0633855994427287e20 6.832356729909674e21 4.381265466255972e23;
                                    0.0 0.0 0.0 0.0 0.0 0.0 -9.150756319887385e8 -2.3425936178911707e11 -3.998026441200931e13 -5.848541536842505e15 -7.985208711635634e17 -1.0550778994470824e20 -1.371936214646111e22 -1.7699057433670775e24 -2.274363584260878e26;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.3610392369296838e11 1.2088520893079981e14 4.126215131504634e16 1.2072126556173556e19 3.296495358272459e21 8.711254830634834e23 2.2654793515098592e26 5.845293495297249e28;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.2135926857366569e14 -1.2427189101943366e17 -8.483627760260006e19 -4.964134186574997e22 -2.7110791504281586e25 -1.4328490580843533e28 -7.452634592906045e30;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.2451508454393354e17 2.5500689314597588e20 3.481694114419724e23 4.074576883618054e26 4.4505245107465464e29 4.704347973165898e32;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.5525616674631408e20 -1.0455292589929025e24 -2.854991896556619e27 -6.682312461883378e30 -1.4597734316732968e34;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0460400207217706e24 8.569159849752745e27 4.679903832611632e31 2.1907298398145424e35;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -8.571252440680644e27 -1.4043139998811168e32 -1.5338853716034814e36;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.4044854458583944e32 4.602217908988787e36;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -4.602498823223603e36]
        extrapolation_scalars = T[-1.0, 0.5, -0.125, 0.015625, -0.0009765625,
            3.0517578125e-5, -4.76837158203125e-7,
            3.725290298461914e-9, -1.4551915228366852e-11,
            2.842170943040401e-14, -2.7755575615628914e-17,
            1.3552527156068805e-20, -3.308722450212111e-24,
            4.0389678347315804e-28, -2.465190328815662e-32,
            7.52316384526264e-37]
        extrapolation_scalars_2 = T[-0.5, 0.125, -0.015625, 0.0009765625, -3.0517578125e-5,
            4.76837158203125e-7, -3.725290298461914e-9,
            1.4551915228366852e-11, -2.842170943040401e-14,
            2.7755575615628914e-17, -1.3552527156068805e-20,
            3.308722450212111e-24, -4.0389678347315804e-28,
            2.465190328815662e-32, -7.52316384526264e-37]
    else # sequence == :bulirsch
        subdividing_sequence = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
        extrapolation_weights = T[-1.0 -2.0 -3.0 -4.0 -4.8 -5.485714285714286 -5.984415584415585 -6.383376623376623 -6.660914737436476 -6.875782954773137 -7.022076209130012 -7.13353773625906 -7.20862760716705 -7.265388454467578 -7.3034271374752615 -7.332068028210459;
                                  0.0 4.0 24.0 96.0 288.0 768.0 1843.2 4213.028571428571 9192.062337662337 19609.732987012987 40924.660146809714 84489.62094825231 172574.54491557917 350627.6468126053 708636.9282949497 1.4284334932559617e6;
                                  0.0 0.0 -27.0 -324.0 -1944.0 -9331.2 -37324.8 -137814.64615384614 -472507.35824175825 -1.5641622893520272e6 -5.005319325926487e6 -1.5754447714391567e7 -4.878796711553518e7 -1.4987663497892407e8 -4.567668875548162e8 -1.3865492871229203e9;
                                  0.0 0.0 0.0 256.0 3072.0 24576.0 147456.0 786432.0 3.7748736e6 1.7256565028571427e7 7.530137467012987e7 3.212858652592208e8 1.3410192636906607e9 5.5371117984646635e9 2.2619690751174793e10 9.191493384604361e10;
                                  0.0 0.0 0.0 0.0 -1555.2 -37324.8 -447897.6 -4.29981696e6 -3.439853568e7 -2.5401995579076922e8 -1.7418511254224176e9 -1.1532255726934628e10 -7.38064366523816e10 -4.6461756843466455e11 -2.877631391595342e12 -1.768016726996178e13;
                                  0.0 0.0 0.0 0.0 0.0 22469.485714285714 539267.6571428571 8.628282514285713e6 1.0353939017142858e8 1.1044201618285713e9 1.0602433553554285e10 9.693653534678204e10 8.459915812082797e11 7.219128159643986e12 6.026402637615849e13 4.976642178160185e14;
                                  0.0 0.0 0.0 0.0 0.0 0.0 -217162.47272727272 -1.0423798690909091e7 -2.501711685818182e8 -4.803286436770909e9 -7.685258298833455e10 -1.1350535333661716e12 -1.5566448457593213e13 -2.061212485419239e14 -2.638351981336626e15 -3.3217283961746376e16;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.663693137582418e6 2.7185727060395604e8 8.699432659326593e9 2.0878638382383826e11 4.454109521575216e12 8.551890281424414e13 1.5637742228890358e15 2.729496825406317e16 4.658341248693448e17;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -9.94469241567476e7 -9.54690471904777e9 -4.5825142651429297e11 -1.7596854778148848e13 -5.630993529007631e14 -1.6633088577991774e16 -4.5622185813920294e17 -1.208201334658303e19;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.95451559685786e9 4.7563349729835455e11 3.044054382709469e13 1.461146103700545e15 6.234223375788993e16 2.393941776302973e18 8.754987067622302e19;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.664005179966794e11 -3.194889945536245e13 -3.067094347714795e15 -2.3555284590449626e17 -1.507538213788776e19 -8.906071909152154e20;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.6222283049151686e13 3.1146783454371235e15 3.986788282159518e17 3.827316750873137e19 3.265976960745077e21;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.066453178967725e15 -4.095180207236064e17 -7.862745997893244e19 -1.2077177852764021e22;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0573083217291824e17 7.90006395544006e19 2.0224163725926554e22;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.676457043743931e19 -2.0555190095953393e22;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0272113415037764e22]
        extrapolation_weights_2 = T[-2.0 -12.0 -48.0 -144.0 -384.0 -921.6 -2106.5142857142855 -4596.0311688311685 -9804.866493506493 -20462.330073404857 -42244.810474126156 -86287.27245778959 -175313.82340630266 -354318.46414747485 -714216.7466279808;
                                    0.0 18.0 216.0 1296.0 6220.8 24883.2 91876.43076923076 315004.9054945055 1.0427748595680182e6 3.336879550617658e6 1.0502965142927712e7 3.2525311410356782e7 9.991775665261604e7 3.0451125836987746e8 9.243661914152801e8;
                                    0.0 0.0 -192.0 -2304.0 -18432.0 -110592.0 -589824.0 -2.8311552e6 -1.2942423771428572e7 -5.64760310025974e7 -2.409643989444156e8 -1.0057644477679955e9 -4.1528338488484974e9 -1.6964768063381096e10 -6.893620038453271e10;
                                    0.0 0.0 0.0 1296.0 31104.0 373248.0 3.5831808e6 2.86654464e7 2.116832964923077e8 1.4515426045186813e9 9.610213105778856e9 6.150536387698467e10 3.8718130702888715e11 2.3980261596627847e12 1.473347272496815e13;
                                    0.0 0.0 0.0 0.0 -19660.8 -471859.2 -7.5497472e6 -9.05969664e7 -9.663676416e8 -9.27712935936e9 -8.481946842843428e10 -7.402426335572446e11 -6.316737139688488e12 -5.273102307913868e13 -4.354561905890162e14;
                                    0.0 0.0 0.0 0.0 0.0 199065.6 9.5551488e6 2.293235712e8 4.40301256704e9 7.044820107264e10 1.0404657389189907e12 1.4269244419460445e13 1.8894447783009694e14 2.4184893162252405e15 3.0449176964934176e16;
                                    0.0 0.0 0.0 0.0 0.0 0.0 -5.309712316483516e6 -2.5486619119120878e8 -8.155718118118681e9 -1.9573723483484836e11 -4.1757276764767646e12 -8.017397138835389e13 -1.466038333958471e15 -2.558903273818422e16 -4.367194920650107e17;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 9.530330231688312e7 9.149117022420778e9 4.391576170761974e11 1.686365249572598e13 5.396368798632314e14 1.594004322057545e16 4.372126140500695e17 1.1578596123808737e19;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -4.799686984456052e9 -4.6076995050778094e11 -2.948927683249798e13 -1.415485287959903e15 -6.0394038952955864e16 -2.3191310957935053e18 -8.481393721759105e19;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.6293384053841528e11 3.128329738337573e13 3.00319654880407e15 2.306454949481526e17 1.4761311676681767e19 8.72052874437815e20;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.596880987650869e13 -3.0660114962896685e15 -3.9244947152507757e17 -3.767514926640745e19 -3.2149460707334354e21;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0553442916868112e15 4.052522080077355e17 7.780842393748523e19 1.195137391679773e22;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.0412356004656733e17 -7.838344705788186e19 -2.0066162446817756e22;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.662517163307765e19 2.0448131814203635e22;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0231987972010274e22]
        extrapolation_scalars = T[-1.0, 0.5, -0.16666666666666666, 0.041666666666666664,
            -0.006944444444444444, 0.0008680555555555555,
            -7.233796296296296e-5, 4.521122685185185e-6,
            -1.8838011188271604e-7, 5.886878496334876e-9,
            -1.226433020069766e-10, 1.9163015938590095e-12,
            -1.9961474936031345e-14, 1.5594902293774489e-16,
            -8.122344944674213e-19, 3.1727909940133645e-21]
        extrapolation_scalars_2 = T[-0.5, 0.16666666666666666, -0.041666666666666664,
            0.006944444444444444, -0.0008680555555555555,
            7.233796296296296e-5, -4.521122685185185e-6,
            1.8838011188271604e-7, -5.886878496334876e-9,
            1.226433020069766e-10, -1.9163015938590095e-12,
            1.9961474936031345e-14, -1.5594902293774489e-16,
            8.122344944674213e-19, -3.1727909940133645e-21]
    end
    extrapolation_coefficients(subdividing_sequence,
        extrapolation_weights, extrapolation_scalars,
        extrapolation_weights_2, extrapolation_scalars_2)
end

function create_extrapolation_coefficients(T::Type{<:CompiledFloats},
        alg::Union{ExtrapolationMidpointDeuflhard,
            ExtrapolationMidpointHairerWanner,
            ImplicitDeuflhardExtrapolation,
            ImplicitHairerWannerExtrapolation})
    # Compute and return extrapolation_coefficients

    @unpack min_order, init_order, max_order, sequence = alg

    max_order > 15 &&
        error("max_order > 15 not allowed for Float32 or Float64 with this algorithm. That's a bad idea.")

    # Initialize subdividing_sequence:
    if sequence == :harmonic
        subdividing_sequence = [
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21
        ]
        extrapolation_weights = T[-1.0 -1.3333333333333333 -1.5 -1.6 -1.6666666666666667 -1.7142857142857142 -1.75 -1.7777777777777777 -1.8 -1.8181818181818181 -1.8333333333333333 -1.8461538461538463 -1.8571428571428572 -1.8666666666666667 -1.875 -1.8823529411764706;
                                  0.0 5.333333333333333 38.4 204.8 975.2380952380952 4388.571428571428 19114.666666666668 81555.91111111111 343170.32727272727 1.4298763636363635e6 5.915044102564103e6 2.433618145054945e7 9.970459794285715e7 4.0712710826666665e8 1.6579836988235295e9 6.737203601568627e9;
                                  0.0 0.0 -72.9 -1499.6571428571428 -21088.928571428572 -253067.14285714287 -2.79006525e6 -2.9219592436363637e7 -2.9584837341818184e8 -2.925972923916084e9 -2.8449861733434067e10 -2.7311867264096704e11 -2.596334381793193e12 -2.449162486354648e13 -2.2960898309574828e14 -2.141777720860745e15;
                                  0.0 0.0 0.0 1872.4571428571428 83220.31746031746 2.3967451428571427e6 5.6940854303030305e7 1.214738225131313e9 2.422001138107972e10 4.613335501158042e11 8.50611193356378e12 1.5311001480414803e14 2.7059443139242895e15 4.7143563158147624e16 8.120422362168969e17 1.3858854164768373e19;
                                  0.0 0.0 0.0 0.0 -77504.96031746031 -6.341314935064935e6 -3.236712831439394e8 -1.3278821872571873e10 -4.80171683784965e11 -1.6005722792832168e13 -5.043469942533053e14 -1.525755612867142e16 -4.4766093502525523e17 -1.2827711003647664e19 -3.607793719775906e20 -9.995618963881296e21;
                                  0.0 0.0 0.0 0.0 0.0 4.711650077922078e6 6.393346721118882e8 5.260811016234965e10 3.4090055385202573e12 1.9175656154176447e14 9.826959789128542e15 4.7169406987817e17 2.15773437679608e19 9.515608601670712e20 4.078117972144591e22 1.708360692331116e24;
                                  0.0 0.0 0.0 0.0 0.0 0.0 -3.952348909376457e8 -8.263044119869713e10 -1.0248756909925902e13 -9.846844874242534e14 -8.108603230469998e16 -6.022558357283822e18 -4.1560671463889437e20 -2.715297202307443e22 -1.7009177076954297e24 -1.0307396968759164e26;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.3741255122091064e10 1.3338509797230594e13 2.371290630618772e15 3.221627130440662e17 3.711314454267642e19 3.8230073464151254e21 3.6330154661690406e23 3.249405137443117e25 2.772825717284793e27;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -6.174193142616171e12 -2.6321560239574205e15 -6.449440297701669e17 -1.1940678036887662e20 -1.8574538823517637e22 -2.5642554640188345e24 -3.245385821648838e26 -3.845504022726303e28;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.082508822446903e15 6.237312738860727e17 2.041302350899874e20 4.99971155510259e22 1.0207744425001122e25 1.837393996500202e27 3.015210660923408e29;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.30788366240374e17 -1.748372388422729e20 -7.448430618928414e22 -2.355293074113417e25 -6.165659032955555e27 -1.4147218829987503e30;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.87960511179021e19 5.723442800021062e22 3.1065086459191243e25 1.2426034583676497e28 4.089940525827236e30;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.7640007344435631e22 -2.1641022343595773e25 -1.4694640618129094e28 -7.307458985088932e30;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 6.155890364150897e24 9.361198795139813e27 7.828458512415587e30;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.472332709198601e27 -4.593753679027078e30;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.1322357960170301e30]
        extrapolation_weights_2 = T[-4.0 -28.8 -153.6 -731.4285714285714 -3291.4285714285716 -14336.0 -61166.933333333334 -257377.74545454545 -1.0724072727272727e6 -4.436283076923077e6 -1.8252136087912086e7 -7.477844845714286e7 -3.053453312e8 -1.2434877741176472e9 -5.052902701176471e9;
                                    0.0 64.8 1333.0285714285715 18745.714285714286 224948.57142857142 2.480058e6 2.5972971054545455e7 2.6297633192727274e8 2.6008648212587414e9 2.5288765985274727e10 2.4277215345863736e11 2.3078527838161714e12 2.1770333212041316e13 2.0409687386288734e14 1.9038024185428845e15;
                                    0.0 0.0 -1755.4285714285713 -78019.04761904762 -2.2469485714285714e6 -5.338205090909091e7 -1.138817086060606e9 -2.2706260669762238e10 -4.325002032335664e11 -7.974479937716044e12 -1.4354063887888878e14 -2.5368227943040215e15 -4.41970904607634e16 -7.612895964533408e17 -1.299267577947035e19;
                                    0.0 0.0 0.0 74404.76190476191 6.087662337662337e6 3.107244318181818e8 1.2747668997668997e10 4.609648164335664e11 1.536549388111888e13 4.8417311448317306e14 1.4647253883524564e16 4.29754497624245e17 1.2314602563501758e19 3.4634819709848696e20 9.595794205326044e21;
                                    0.0 0.0 0.0 0.0 -4.580770909090909e6 -6.215753756643356e8 -5.114677376895105e10 -3.314310940228028e12 -1.8642999038782656e14 -9.553988683874972e15 -4.585914568259986e17 -2.0977973107739664e19 -9.251286140513193e20 -3.964836917362797e22 -1.6609062286552515e24;
                                    0.0 0.0 0.0 0.0 0.0 3.8716887275524473e8 8.094410566402983e10 1.0039598605641701e13 9.64588885640085e14 7.943121531888978e16 5.899649003053539e18 4.071249449523863e20 2.659882973688924e22 1.6662051014159312e24 1.0097041928580407e26;
                                    0.0 0.0 0.0 0.0 0.0 0.0 -4.3057798010808395e10 -1.3130095581648865e13 -2.3342392145153535e15 -3.1712892065275264e17 -3.65332516591971e19 -3.763272856627389e21 -3.5762495995101494e23 -3.198633182170568e25 -2.729500315452218e27;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 6.097968535917206e12 2.59966027057523e15 6.369817577976957e17 1.1793262258654482e20 1.8345223529400134e22 2.532597989154405e24 3.205319330023544e26 3.7980286644210397e28;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0716837342224339e15 -6.174939611472119e17 -2.0208893273908753e20 -4.949714439551565e22 -1.010566698075111e25 -1.8190200565352e27 -2.9850585543141743e29;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.2888102437061885e17 1.733923029840723e20 7.38687334108603e22 2.3358278420959505e25 6.114703173179063e27 1.4030299666103307e30;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -5.838774520736112e19 -5.6836966694653606e22 -3.0849356692113527e25 -1.233974267684541e28 -4.061538161064547e30;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.753562860275258e22 2.1512968956947278e25 1.4607690081927147e28 7.2642195828103e30;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -6.124482760252167e24 -9.313437576797262e27 -7.788517397556324e30;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.461344563824385e27 4.57333699600918e30;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.1278129999388386e30]
        extrapolation_scalars = T[-1.0, 0.25, -0.027777777777777776, 0.001736111111111111,
            -6.944444444444444e-5, 1.9290123456790124e-6,
            -3.936759889140842e-8, 6.151187326782565e-10,
            -7.594058428126624e-12, 7.594058428126623e-14,
            -6.276081345559193e-16, 4.358389823304995e-18,
            -2.5789288895295828e-20, 1.3157800456783586e-22,
            -5.8479113141260385e-25, 2.2843403570804838e-27]
        extrapolation_scalars_2 = T[-0.25, 0.027777777777777776, -0.001736111111111111,
            6.944444444444444e-5, -1.9290123456790124e-6,
            3.936759889140842e-8, -6.151187326782565e-10,
            7.594058428126624e-12, -7.594058428126623e-14,
            6.276081345559193e-16, -4.358389823304995e-18,
            2.5789288895295828e-20, -1.3157800456783586e-22,
            5.8479113141260385e-25, -2.2843403570804838e-27]
    elseif sequence == :romberg
        subdividing_sequence = [
            1,
            2,
            4,
            8,
            16,
            32,
            64,
            128,
            256,
            512,
            1024,
            2048,
            4096,
            8192,
            16384,
            32768
        ]
        extrapolation_weights = T[-1.0 -1.3333333333333333 -1.4222222222222223 -1.4447971781305116 -1.4504630494172979 -1.451880901860521 -1.452235451531305 -1.4523240943593299 -1.4523462554044868 -1.4523517956869105 -1.4523531807588372 -1.4523535270269017 -1.4523536135939228 -1.4523536352356785 -1.4523536406461173 -1.452353641998727;
                                  0.0 5.333333333333333 28.444444444444443 121.36296296296297 493.15743680188126 1980.3655501377505 7929.205565360925 31724.5675171852 126906.01579724405 507631.8090356717 2.0305349820189304e6 8.122147673959353e6 3.248859844172289e7 1.2995440151277749e8 5.19817613796996e8 2.07927046293387e9;
                                  0.0 0.0 -91.02222222222223 -1941.8074074074075 -33140.17975308642 -538659.4296374682 -8.652349112921841e6 -1.3857291091506496e8 -2.2177080072599993e9 -3.54854939788296e10 -5.677765672441477e11 -9.084459730370057e12 -1.4535149430390788e14 -2.325624463334606e15 -3.720999363124215e16 -5.953599069714284e17;
                                  0.0 0.0 0.0 5917.889241622575 504993.2152851264 3.4474203496797964e7 2.241370436871182e9 1.4401024799097037e11 9.225665310201598e12 5.905867660750885e14 3.77998601491761e16 2.4192279640364677e18 1.5483118033241417e20 9.909204991428803e21 6.341892706539483e23 4.05881157410929e25;
                                  0.0 0.0 0.0 0.0 -1.5209207425057925e6 -5.1914094677531046e8 -1.4176008786611145e11 -3.6866623485688414e13 -9.474866810815984e15 -2.4279369357326935e18 -6.217036386624773e20 -1.591658462098873e23 -4.074707838080507e25 -1.0431291857706623e28 -2.670413262277438e30 -6.836259581321538e32;
                                  0.0 0.0 0.0 0.0 0.0 1.5589452477944808e9 2.1284799116553977e12 2.3248676581708025e15 2.4184528070774876e18 2.486207422190278e21 2.5483650380553205e24 2.6101630458060033e27 2.6729701040532997e30 2.7371631523457504e33 2.8028657600863993e36 2.870137275504677e39;
                                  0.0 0.0 0.0 0.0 0.0 0.0 -6.386999060908798e12 -3.4881530871309916e16 -1.5239973381214444e20 -6.341377114357268e23 -2.607614058456583e27 -1.069122783562139e31 -4.380196305334128e34 -1.7942379182565492e38 -7.349310654759861e41 -3.010289127531343e45;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0465098000284594e17 2.2861355418221703e21 3.995311436502874e25 6.649821721972122e29 1.0937793665786103e34 1.7937998718897874e38 2.93967940527153e42 4.816664723480877e46 7.891743901406595e50;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -6.858511278043386e21 -5.993058601571351e26 -4.189451610800854e31 -2.7891799442838834e36 -1.835085270122301e41 -1.2038170852496654e46 -7.891262227584487e50 -5.171933283225826e55;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.7979244390088466e27 6.284201388527134e32 1.7571900680469942e38 4.679485289631606e43 1.2315095760466199e49 3.231484209329825e54 8.473210620642257e59;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.8852622144842938e33 -2.635787615753444e39 -2.948078543974702e45 -3.140352413792322e51 -3.305811498811932e57 -3.46978305819599e63;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 7.907364732522996e39 4.422118870277351e46 1.9784224923818425e53 8.429821331796452e59 3.549588914822444e66;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.3266357401568573e47 -2.9676339154575293e54 -5.310787755579382e61 -9.051452272305827e68;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 8.902901879036164e54 7.966181752074431e62 5.702415016525277e70;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.389854534525231e63 -8.553622556652643e71;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.5660867693856473e72]
        extrapolation_weights_2 = T[-4.0 -21.333333333333332 -91.02222222222223 -369.86807760141096 -1485.274162603313 -5946.904174020694 -23793.4256378889 -95179.51184793304 -380723.8567767538 -1.522901236514198e6 -6.091610755469514e6 -2.4366448831292167e7 -9.746580113458312e7 -3.89863210347747e8 -1.5594528472004025e9;
                                    0.0 85.33333333333333 1820.4444444444443 31068.91851851852 504993.2152851264 8.111577293364226e6 1.299121039828734e8 2.0791012568062494e9 3.3267650605152744e10 5.322905317913885e11 8.516680997221928e12 1.3626702590991364e14 2.1802729343761932e15 3.4884369029289516e16 5.5814991278571405e17;
                                    0.0 0.0 -5825.422222222222 -497102.6962962963 -3.3935544067160495e7 -2.2063490237950697e9 -1.4176008786611145e11 -9.081514289729697e12 -5.813588478551652e14 -3.7209237334345224e16 -2.3814275270983977e18 -1.524119431397202e20 -9.754373663437728e21 -6.242800632999802e23 -3.995392643263833e25;
                                    0.0 0.0 0.0 1.5149796458553793e6 5.1711305245196944e8 1.4120633752288446e11 3.6722613237697445e13 9.437855612336234e15 2.4184528070774876e18 6.19275108823952e20 1.585441046231299e23 4.058791010588005e25 1.0390544623887458e28 2.659981960471667e30 6.809555442332001e32;
                                    0.0 0.0 0.0 0.0 -1.5574228403259315e9 -2.1264013179916716e12 -2.32259727959837e15 -2.416091036758076e18 -2.4837794852545453e21 -2.545876400322845e24 -2.607614058456583e27 -2.6703597816860604e30 -2.7344901414547876e33 -2.8001285864925646e36 -2.867334407071567e39;
                                    0.0 0.0 0.0 0.0 0.0 6.385439734966193e12 3.4873014872562036e16 1.523625268458817e20 6.339828926585209e23 2.606977433930593e27 1.0688617672575583e31 4.379126921470521e34 1.7937998718897874e38 7.34751638946329e41 3.009554193662317e45;
                                    0.0 0.0 0.0 0.0 0.0 0.0 -1.0464459261392974e17 -2.2859960071821667e21 -3.995067582045079e25 -6.649415849064287e29 -1.093712607584068e34 -1.7936903870343255e38 -2.9394999814797044e42 -4.816370737596875e46 -7.891262227584487e50;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 6.858406625466511e21 5.99296715475431e26 4.1893876848424375e31 2.789137384775456e36 1.8350572689432526e41 1.2037987164586916e46 7.89114181647872e50 5.171854365786813e55;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.7979175804714053e27 -6.284177416201281e32 -1.7571833648988464e38 -4.679467438811868e43 -1.2315048782104076e49 -3.231471882195848e54 -8.47317829790887e59;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.8852604165581403e33 2.6357851020704912e39 2.948075732467912e45 3.140349418918881e51 3.305808346144411e57 3.4697797491530044e63;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -7.907362847260331e39 -4.4221178159620535e46 -1.978422020689163e53 -8.429819321970427e59 -3.549588068534498e66;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.3266356610832054e47 2.967633738572764e54 5.310787439031764e61 9.051451732797232e68;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -8.902901746372587e54 -7.966181633369073e62 -5.702414931552672e70;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.3898545256223295e63 8.553622524787915e71;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.5660867669957927e72]
        extrapolation_scalars = T[-1.0, 0.25, -0.015625, 0.000244140625,
            -9.5367431640625e-7, 9.313225746154785e-10,
            -2.2737367544323206e-13, 1.3877787807814457e-17,
            -2.117582368135751e-22, 8.077935669463161e-28,
            -7.703719777548943e-34, 1.8367099231598242e-40,
            -1.0947644252537633e-47, 1.6313261169996311e-55,
            -6.077163357286271e-64, 5.659799424266695e-73]
        extrapolation_scalars_2 = T[-0.25, 0.015625, -0.000244140625, 9.5367431640625e-7,
            -9.313225746154785e-10, 2.2737367544323206e-13,
            -1.3877787807814457e-17, 2.117582368135751e-22,
            -8.077935669463161e-28, 7.703719777548943e-34,
            -1.8367099231598242e-40, 1.0947644252537633e-47,
            -1.6313261169996311e-55, 6.077163357286271e-64,
            -5.659799424266695e-73]
    else # sequence == :bulirsch
        subdividing_sequence = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
        extrapolation_weights = T[-1.0 -1.3333333333333333 -1.5 -1.6 -1.6457142857142857 -1.6718367346938776 -1.6835279006707577 -1.6901299708694666 -1.693069327340544 -1.6947243315705933 -1.6954602083971546 -1.695874240194077 -1.6960582742950205 -1.6961617997954963 -1.6962078123772122 -1.6962336948493626;
                                  0.0 5.333333333333333 38.4 204.8 921.6 3932.16 16178.029714285714 65739.29534693877 264796.0427960611 1.063337834600653e6 4.260748471165052e6 1.7059653702729277e7 6.82682451256418e7 2.7313966499109036e8 1.0926772230310967e9 4.3709756753076935e9;
                                  0.0 0.0 -72.9 -1499.6571428571428 -17995.885714285716 -188466.0031168831 -1.809273629922078e6 -1.687678722000189e7 -1.5430205458287445e8 -1.4010322512667694e9 -1.2658738458504457e10 -1.1417952888042778e11 -1.0286202719081353e12 -9.262670584090748e12 -8.338439277458397e13 -7.505626090600242e14;
                                  0.0 0.0 0.0 1872.4571428571428 53926.76571428571 1.1504376685714286e6 2.0707878034285713e7 3.5341445178514284e8 5.816192120806923e9 9.4536202090576e10 1.5231567106062036e12 2.4466077986835336e13 3.9213804300291206e14 6.280341834369219e15 1.0052910177255184e17 1.608858416060063e18;
                                  0.0 0.0 0.0 0.0 -57586.834285714285 -4.738573792653061e6 -2.2745154204734695e8 -9.528151870492496e9 -3.6588103182691187e11 -1.36516582563434e13 -4.992606448034158e14 -1.8132753113333124e16 -6.553390301665806e17 -2.3644157580681015e19 -8.52021725370699e20 -3.0689640436338757e22;
                                  0.0 0.0 0.0 0.0 0.0 5.09977563903525e6 5.874941536168609e8 5.013283444197213e10 3.609564079821993e12 2.464129078491814e14 1.6221009705271826e16 1.0546231071871968e18 6.7967878012847596e19 4.36700279749998e21 2.7997424543832918e23 1.7935867203368857e25;
                                  0.0 0.0 0.0 0.0 0.0 0.0 -5.70060368320064e8 -1.8763129837277533e11 -3.602520928757287e13 -6.036515068986755e15 -9.272087145963656e17 -1.3838308524243085e20 -2.0243468469749884e22 -2.9409072775127474e24 -4.251513956009017e26 -6.1356617753586065e28;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.9561237425512534e11 9.013818205676177e13 3.0767166142041348e16 8.860943848907908e18 2.419628400341786e21 6.371227239299971e23 1.6569236045823926e26 4.271386836316456e28 1.0977631674699422e31;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -8.554159840110055e13 -1.1262162440922038e17 -8.649340754628125e19 -5.797259955974749e22 -3.561836516950886e25 -2.126373139447408e28 -1.2442320541680835e31 -7.230324405099858e33;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.1651370210969312e17 2.1475805572858636e20 2.9321633208809658e23 3.377852145654872e26 3.689515303627295e29 3.8860083472261895e32 4.0424356038700874e35;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -2.0269737017532947e20 -1.0674622648776207e24 -3.279244077704051e27 -8.791712994944155e30 -2.1606513856374756e34 -5.159530537984771e37;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.1022516982215851e24 8.126681320648103e27 4.438251558583284e31 2.0451463181951773e35 8.935380607282608e38;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -7.659880204931148e27 -1.6135647078547534e32 -1.9827483130119208e36 -2.126313710861715e40;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.665360734159382e32 4.911348648324117e36 1.0729004833885644e41;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -4.62766803500927e36 -3.8992995301161535e41;
                                  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.023990816545532e41]
        extrapolation_weights_2 = T[-4.0 -28.8 -153.6 -691.2 -2949.12 -12133.522285714285 -49304.47151020408 -198597.0320970458 -797503.3759504899 -3.1955613533737888e6 -1.279474027704696e7 -5.120118384423134e7 -2.0485474874331778e8 -8.195079172733225e8 -3.27823175648077e9;
                                    0.0 64.8 1333.0285714285715 15996.342857142858 167525.33610389612 1.6082432265974027e6 1.5001588640001683e7 1.3715738185144395e8 1.2453620011260173e9 1.1252211963115074e10 1.0149291456038025e11 9.143291305850092e11 8.233484963636221e12 7.411946024407464e13 6.671667636089105e14;
                                    0.0 0.0 -1755.4285714285713 -50556.34285714286 -1.0785353142857142e6 -1.941363565714286e7 -3.313260485485714e8 -5.45268011325649e9 -8.862768945991501e10 -1.427959416193316e12 -2.2936948112658125e13 -3.6762941531523006e14 -5.887820469721143e15 -9.424603291176736e16 -1.508304765056309e18;
                                    0.0 0.0 0.0 55987.2 4.606946742857143e6 2.2113344365714285e8 9.263480985201038e9 3.557176698317199e11 1.327244552700053e13 4.853922935588765e14 1.7629065526851648e16 6.37135168217509e17 2.29873754256621e19 8.283544552215128e20 2.9837150424218233e22;
                                    0.0 0.0 0.0 0.0 -5.020091644675325e6 -5.783145574665974e8 -4.9349508903816315e10 -3.5531646410747744e12 -2.4256270616403794e14 -1.5967556428626954e16 -1.0381446211373969e18 -6.690587991889685e19 -4.2987683787890434e21 -2.755996478533553e23 -1.765561927831622e25;
                                    0.0 0.0 0.0 0.0 0.0 5.661016157622857e8 1.863283032451866e11 3.577503422307583e13 5.994594825452124e15 9.207697651894463e17 1.3742209159491396e20 2.0102888827598844e22 2.920484310307798e24 4.221989553536732e26 6.093053013029727e28;
                                    0.0 0.0 0.0 0.0 0.0 0.0 -1.9484826341819125e11 -8.978607978310253e13 -3.0646981899298996e16 -8.826330786998111e18 -2.410176726902951e21 -6.346339632896457e23 -1.6504512467519928e26 -4.254701731487095e28 -1.0934750300970127e31;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 8.539308868165419e13 1.1242610075573214e17 8.63432453804023e19 5.787195268551182e22 3.555652772997846e25 2.1226815194136454e28 1.2420719290740417e31 7.217771758563227e33;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.1639991919747662e17 -2.145483310647889e20 -2.929299880137918e23 -3.3745534619188814e26 -3.685912261338597e29 -3.882213417199601e32 -4.0384879128506835e35;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0260939388619086e20 1.0669989566029343e24 3.277820794684214e27 8.787897147290099e30 2.1597136029180146e34 5.157291158410993e37;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.1019825938030739e24 -8.124697267591304e27 -4.437168001073864e31 -2.0446470148948364e35 -8.933199117876534e38;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 7.6590490547353e27 1.6133896248786405e32 1.9825331710508737e36 2.126082991058019e40;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.66525908860676e32 -4.911048883391968e36 -1.07283499873992e41;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.627542501479674e36 3.8991937548467816e41;
                                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -4.023929415318473e41]
        extrapolation_scalars = T[-1.0, 0.25, -0.027777777777777776, 0.001736111111111111,
            -4.8225308641975306e-5, 7.535204475308642e-7,
            -5.232780885631001e-9, 2.0440550334496098e-11,
            -3.548706655294462e-14, 3.465533843060998e-17,
            -1.5041379527174468e-20, 3.672211798626579e-24,
            -3.9846048162180767e-28, 2.4320097755237285e-32,
            -6.597248740027475e-37, 1.0066602691692314e-41]
        extrapolation_scalars_2 = T[-0.25, 0.027777777777777776, -0.001736111111111111,
            4.8225308641975306e-5, -7.535204475308642e-7,
            5.232780885631001e-9, -2.0440550334496098e-11,
            3.548706655294462e-14, -3.465533843060998e-17,
            1.5041379527174468e-20, -3.672211798626579e-24,
            3.9846048162180767e-28, -2.4320097755237285e-32,
            6.597248740027475e-37, -1.0066602691692314e-41]
    end
    extrapolation_coefficients(subdividing_sequence,
        extrapolation_weights, extrapolation_scalars,
        extrapolation_weights_2, extrapolation_scalars_2)
end

function generate_sequence(T, alg::ImplicitEulerExtrapolation)
    # Compute and return extrapolation_coefficients

    @unpack min_order, init_order, max_order, sequence = alg

    # Initialize subdividing_sequence:
    if sequence == :harmonic
        subdividing_sequence = BigInt.(1:(max_order + 1))
    elseif sequence == :romberg
        subdividing_sequence = BigInt(2) .^ (0:max_order)
    else # sequence == :bulirsch
        subdividing_sequence = [n == 0 ? BigInt(1) :
                                (isodd(n) ? BigInt(2)^((n + 1) ÷ 2) :
                                 3 * BigInt(2)^(n ÷ 2 - 1)) for n in 0:max_order]
    end

    subdividing_sequence
end

function generate_sequence(T::Type{<:CompiledFloats}, alg::ImplicitEulerExtrapolation)
    # Compute and return extrapolation_coefficients

    @unpack min_order, init_order, max_order, sequence = alg

    # Initialize subdividing_sequence:
    if sequence == :harmonic
        subdividing_sequence = Int.(1:(max_order + 1))
    elseif sequence == :romberg
        subdividing_sequence = Int(2) .^ (0:max_order)
    else # sequence == :bulirsch
        subdividing_sequence = [n == 0 ? Int(1) :
                                (isodd(n) ? Int(2)^((n + 1) ÷ 2) : 3 * Int(2)^(n ÷ 2 - 1))
                                for n in 0:max_order]
    end

    subdividing_sequence
end

@cache mutable struct ExtrapolationMidpointDeuflhardConstantCache{QType,
    extrapolation_coefficients
} <:
                      OrdinaryDiffEqConstantCache
    # Values that are mutated
    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n + alg.min_order - 1)
    n_curr::Int # Storage for the current extrapolation order
    n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter

    # Constant values
    coefficients::extrapolation_coefficients
    stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n + alg.min_order - 1)
end

function alg_cache(alg::ExtrapolationMidpointDeuflhard, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # Initialize cache's members
    QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl

    Q = fill(zero(QType), alg.max_order - alg.min_order + 1)
    n_curr = alg.init_order
    n_old = alg.init_order
    sequence_factor = alg.sequence_factor

    coefficients = create_extrapolation_coefficients(constvalue(uBottomEltypeNoUnits), alg)
    stage_number = Vector{Int}(undef, alg.max_order - alg.min_order + 1)
    for n in 1:length(stage_number)
        s = zero(eltype(coefficients.subdividing_sequence))
        for i in 1:(alg.min_order + n)
            s += coefficients.subdividing_sequence[i]
        end
        stage_number[n] = sequence_factor * Int(s) - alg.min_order - n + 3 - sequence_factor
    end

    # Initialize cache
    ExtrapolationMidpointDeuflhardConstantCache(Q, n_curr, n_old, coefficients,
        stage_number)
end

@cache mutable struct ExtrapolationMidpointDeuflhardCache{uType, uNoUnitsType, rateType,
    QType, extrapolation_coefficients
} <: ExtrapolationMutableCache
    # Values that are mutated
    utilde::uType
    u_temp1::uType
    u_temp2::uType
    u_temp3::Array{uType, 1}
    u_temp4::Array{uType, 1}
    tmp::uType # for get_tmp_cache()
    T::Array{uType, 1}  # Storage for the internal discretisations obtained by the explicit midpoint rule
    res::uNoUnitsType # Storage for the scaled residual of u and utilde

    fsalfirst::rateType
    k::rateType
    k_tmps::Array{rateType, 1}

    # Constant values
    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n + alg.min_order - 1)
    n_curr::Int # Storage for the current extrapolation order
    n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter
    coefficients::extrapolation_coefficients
    stage_number::Vector{Int} # Stage_number[n] contains information for extrapolation order (n + alg.min_order - 1)
end

function alg_cache(alg::ExtrapolationMidpointDeuflhard, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # Initialize cache's members
    utilde = zero(u)
    u_temp1 = zero(u)
    u_temp2 = zero(u)
    u_temp3 = Array{typeof(u), 1}(undef, Threads.nthreads())
    u_temp4 = Array{typeof(u), 1}(undef, Threads.nthreads())

    for i in 1:Threads.nthreads()
        u_temp3[i] = zero(u)
        u_temp4[i] = zero(u)
    end

    tmp = zero(u)
    T = Vector{typeof(u)}(undef, alg.max_order + 1)
    for i in 1:(alg.max_order + 1)
        T[i] = zero(u)
    end
    res = uEltypeNoUnits.(zero(u))

    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    k_tmps = Array{typeof(k), 1}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        k_tmps[i] = zero(rate_prototype)
    end

    cc = alg_cache(alg::ExtrapolationMidpointDeuflhard, u, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p,
        calck, Val(false))
    # Initialize cache
    ExtrapolationMidpointDeuflhardCache(utilde, u_temp1, u_temp2, u_temp3, u_temp4, tmp, T,
        res, fsalfirst, k, k_tmps, cc.Q, cc.n_curr,
        cc.n_old, cc.coefficients, cc.stage_number)
end

@cache mutable struct ImplicitDeuflhardExtrapolationConstantCache{QType,
    extrapolation_coefficients,
    TF, UF} <:
                      OrdinaryDiffEqConstantCache
    # Values that are mutated
    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n + alg.min_order - 1)
    n_curr::Int # Storage for the current extrapolation order
    n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter

    # Constant values
    coefficients::extrapolation_coefficients
    stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n + alg.min_order - 1)

    tf::TF
    uf::UF
end

@cache mutable struct ImplicitDeuflhardExtrapolationCache{uType, QType,
    extrapolation_coefficients,
    rateType, JType, WType, F, JCType,
    GCType, uNoUnitsType, TFType,
    UFType} <:
                      ExtrapolationMutableCache
    # Values that are mutated
    utilde::uType
    u_temp1::uType
    u_temp2::uType
    u_temp3::Array{uType, 1}
    u_temp4::Array{uType, 1}
    tmp::uType # for get_tmp_cache()
    T::Array{uType, 1}  # Storage for the internal discretisations obtained by the explicit midpoint rule
    res::uNoUnitsType # Storage for the scaled residual of u and utilde

    fsalfirst::rateType
    k::rateType
    k_tmps::Array{rateType, 1}

    # Constant values
    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n + alg.min_order - 1)
    n_curr::Int # Storage for the current extrapolation order
    n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter
    coefficients::extrapolation_coefficients
    stage_number::Vector{Int} # Stage_number[n] contains information for extrapolation order (n + alg.min_order - 1)

    du1::rateType
    du2::rateType
    J::JType
    W::WType
    tf::TFType
    uf::UFType
    linsolve_tmps::Array{rateType, 1}
    linsolve::Array{F, 1}
    jac_config::JCType
    grad_config::GCType
    # Values to check overflow in T1 computation
    diff1::Array{uType, 1}
    diff2::Array{uType, 1}
end

function alg_cache(alg::ImplicitDeuflhardExtrapolation, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # Initialize cache's members
    QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl

    Q = fill(zero(QType), alg.max_order - alg.min_order + 1)
    n_curr = alg.init_order
    n_old = alg.init_order

    coefficients = create_extrapolation_coefficients(constvalue(uBottomEltypeNoUnits), alg)
    stage_number = Vector{Int}(undef, alg.max_order - alg.min_order + 1)

    #==
    Work calculation in Deuflhard is referenced from here: https://link.springer.com/article/10.1007/BF01418332
    A[1] := CJAC + CLR + (N[1] + 1)(CF + CS)
    A[J] := A[J-1] - N[J]*(CF + CS) + CLR + CS      J = 2, 3, 4.....
    CF = 1; CJ = n*CF ; CS = CLR = 0
    n = Dimension of the jacobian (particularly gaussian decomposition of I - hJ (n,n) matrix)
    Since we are using 4*N sequence and doing 4*N - 1 Computations
    A[J] := A[J-1] - (4*N[J] - 1)*(CF + CS) + CLR + CS      J = 2, 3, 4.....
    ===#
    for n in 1:length(stage_number)
        s = zero(eltype(coefficients.subdividing_sequence))
        for i in 1:(alg.min_order + n)
            s += coefficients.subdividing_sequence[i]
        end
        stage_number[n] = 4 * Int(s) - alg.min_order - n - 1
    end

    #Update stage_number by the jacobian size
    jac_dim = rate_prototype isa Union{CompiledFloats, BigFloat} ? 1 :
              sum(size(rate_prototype))
    stage_number = stage_number .+ jac_dim

    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    ImplicitDeuflhardExtrapolationConstantCache(Q, n_curr, n_old, coefficients,
        stage_number, tf, uf)
end

function alg_cache(alg::ImplicitDeuflhardExtrapolation, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    utilde = zero(u)
    u_temp1 = zero(u)
    u_temp2 = zero(u)
    u_temp3 = Array{typeof(u), 1}(undef, Threads.nthreads())
    u_temp4 = Array{typeof(u), 1}(undef, Threads.nthreads())

    for i in 1:Threads.nthreads()
        u_temp3[i] = zero(u)
        u_temp4[i] = zero(u)
    end

    tmp = zero(u)
    T = Vector{typeof(u)}(undef, alg.max_order + 1)
    for i in 1:(alg.max_order + 1)
        T[i] = zero(u)
    end
    res = uEltypeNoUnits.(zero(u))

    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    k_tmps = Array{typeof(k), 1}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        k_tmps[i] = zero(rate_prototype)
    end

    cc = alg_cache(alg::ImplicitDeuflhardExtrapolation, u, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p,
        calck, Val(false))

    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)

    if DiffEqBase.has_jac(f) && !DiffEqBase.has_Wfact(f) && f.jac_prototype !== nothing
        W_el = WOperator(f, dt, true)
        J = nothing # is J = W.J better?
    else
        J = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
        W_el = zero(J)
    end

    W = Array{typeof(W_el), 1}(undef, Threads.nthreads())
    W[1] = W_el
    for i in 2:Threads.nthreads()
        if W_el isa WOperator
            W[i] = WOperator(f, dt, true)
        else
            W[i] = zero(W_el)
        end
    end
    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)
    linsolve_tmp = zero(rate_prototype)
    linsolve_tmps = Array{typeof(linsolve_tmp), 1}(undef, Threads.nthreads())

    for i in 1:Threads.nthreads()
        linsolve_tmps[i] = zero(rate_prototype)
    end

    linprob = LinearProblem(W[1], _vec(linsolve_tmps[1]); u0 = _vec(k_tmps[1]))
    linsolve1 = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true))
    #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
    #Pr = Diagonal(_vec(weight)))

    linsolve = Array{typeof(linsolve1), 1}(undef, Threads.nthreads())
    linsolve[1] = linsolve1
    for i in 2:Threads.nthreads()
        linprob = LinearProblem(W[i], _vec(linsolve_tmps[i]); u0 = _vec(k_tmps[i]))
        linsolve[i] = init(linprob, alg.linsolve,
            alias = LinearAliasSpecifier(alias_A = true, alias_b = true))
        #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
        #Pr = Diagonal(_vec(weight)))
    end
    grad_config = build_grad_config(alg, f, tf, du1, t)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, du1, du2)

    diff1 = Array{typeof(u), 1}(undef, Threads.nthreads())
    diff2 = Array{typeof(u), 1}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        diff1[i] = zero(u)
        diff2[i] = zero(u)
    end

    ImplicitDeuflhardExtrapolationCache(utilde, u_temp1, u_temp2, u_temp3, u_temp4, tmp, T,
        res, fsalfirst, k, k_tmps, cc.Q, cc.n_curr,
        cc.n_old, cc.coefficients, cc.stage_number,
        du1, du2, J, W, tf, uf, linsolve_tmps, linsolve,
        jac_config, grad_config, diff1, diff2)
end

@cache mutable struct ExtrapolationMidpointHairerWannerConstantCache{QType,
    extrapolation_coefficients
} <:
                      OrdinaryDiffEqConstantCache
    # Values that are mutated
    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
    n_curr::Int # Storage for the current extrapolation order
    n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter

    # Constant values
    coefficients::extrapolation_coefficients
    stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n - 1)
    sigma::Rational{Int} # Parameter for order selection

    #Stepsizing caches
    work::Array{QType, 1}
    dt_new::Array{QType, 1}
end

function alg_cache(alg::ExtrapolationMidpointHairerWanner, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # Initialize cache's members
    QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl

    Q = fill(zero(QType), alg.max_order + 1)
    n_curr = alg.init_order
    n_old = alg.init_order
    sequence_factor = alg.sequence_factor

    coefficients = create_extrapolation_coefficients(constvalue(uBottomEltypeNoUnits), alg)
    stage_number = Vector{Int}(undef, alg.max_order + 1)
    for n in 1:length(stage_number)
        s = zero(eltype(coefficients.subdividing_sequence))
        for i in 1:n
            s += coefficients.subdividing_sequence[i]
        end
        stage_number[n] = sequence_factor * Int(s) - n + 3 - sequence_factor
    end
    sigma = 9 // 10

    work = fill(zero(eltype(Q)), alg.max_order + 1)
    dt_new = fill(zero(eltype(Q)), alg.max_order + 1)
    # Initialize the constant cache
    ExtrapolationMidpointHairerWannerConstantCache(Q, n_curr, n_old, coefficients,
        stage_number, sigma, work, dt_new)
end

@cache mutable struct ExtrapolationMidpointHairerWannerCache{uType, uNoUnitsType, rateType,
    QType,
    extrapolation_coefficients} <:
                      ExtrapolationMutableCache
    # Values that are mutated
    utilde::uType
    u_temp1::uType
    u_temp2::uType
    u_temp3::Array{uType, 1}
    u_temp4::Array{uType, 1}
    tmp::uType # for get_tmp_cache()
    T::Array{uType, 1}  # Storage for the internal discretisations obtained by the explicit midpoint rule
    res::uNoUnitsType # Storage for the scaled residual of u and utilde

    fsalfirst::rateType
    k::rateType
    k_tmps::Array{rateType, 1}

    # Constant values
    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
    n_curr::Int # Storage for the current extrapolation order
    n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter
    coefficients::extrapolation_coefficients
    stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n - 1)
    sigma::Rational{Int} # Parameter for order selection

    #Stepsizing caches
    work::Array{QType, 1}
    dt_new::Array{QType, 1}
end

function alg_cache(alg::ExtrapolationMidpointHairerWanner, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # Initialize cache's members
    utilde = zero(u)
    u_temp1 = zero(u)
    u_temp2 = zero(u)
    u_temp3 = Array{typeof(u), 1}(undef, Threads.nthreads())
    u_temp4 = Array{typeof(u), 1}(undef, Threads.nthreads())

    for i in 1:Threads.nthreads()
        u_temp3[i] = zero(u)
        u_temp4[i] = zero(u)
    end
    tmp = zero(u)
    T = Vector{typeof(u)}(undef, alg.max_order + 1)
    for i in 1:(alg.max_order + 1)
        T[i] = zero(u)
    end
    res = uEltypeNoUnits.(zero(u))
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    k_tmps = Array{typeof(k), 1}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        k_tmps[i] = zero(rate_prototype)
    end

    cc = alg_cache(alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(false))

    # Initialize the cache
    ExtrapolationMidpointHairerWannerCache(utilde, u_temp1, u_temp2, u_temp3, u_temp4, tmp,
        T, res, fsalfirst, k, k_tmps,
        cc.Q, cc.n_curr, cc.n_old, cc.coefficients,
        cc.stage_number, cc.sigma, cc.work, cc.dt_new)
end

@cache mutable struct ImplicitHairerWannerExtrapolationConstantCache{QType,
    extrapolation_coefficients,
    TF, UF} <:
                      OrdinaryDiffEqConstantCache
    # Values that are mutated
    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
    n_curr::Int # Storage for the current extrapolation order
    n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter

    # Constant values
    coefficients::extrapolation_coefficients
    stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n - 1)
    sigma::Rational{Int} # Parameter for order selection

    tf::TF
    uf::UF

    #Stepsizing caches
    work::Array{QType, 1}
    dt_new::Array{QType, 1}
end

function alg_cache(alg::ImplicitHairerWannerExtrapolation, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # Initialize cache's members
    QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl

    Q = fill(zero(QType), alg.max_order + 1)
    n_curr = alg.init_order
    n_old = alg.init_order

    coefficients = create_extrapolation_coefficients(constvalue(uBottomEltypeNoUnits), alg)
    #==Work Calculation (A[J] denotes Jth order work)
    Default values are used from https://github.com/luchr/ODEInterface.jl/blob/master/src/Seulex.jl#L393-L399

    ║ WKFCN      │ estimated works (complexity)        │     1.0 ║
    ║ WKJAC      │ for a call to                       │     5.0 ║
    ║ WKDEC      │ WKFCN: right-hand side f            │     1.0 ║
    ║ WKSOL      │ WKJAC: JACOBIMATRIX                 │     1.0 ║
    ║ WKROW      │ WKDEC: LU-decomposition             │     2.0 ║
    ║            │ WKSOL: Forward- and Backward subst. │         ║
    ║            | WKROW: Tot. work in one iteration   |         ║
    ╚════════════╧═════════════════════════════════════╧═════════╝
    WKROW = WKFCN + WKSOL
    A[1] = WKJAC + (N[1] + 1)* WKROw + WKDEC
    A[J] = A[J - 1] + N[J]* WKROW + WKDEC

    Since we are using 4*N Sequence and only performing 4*N - 1 computations, The modified Work Equation becomes:
    A[J] = A[J - 1] + (4*N[J] - 1)* WKROW + WKDEC
    ==#
    stage_number = Vector{Int}(undef, alg.max_order + 1)
    for n in 1:length(stage_number)
        s = zero(eltype(coefficients.subdividing_sequence))
        for i in 1:n
            s += coefficients.subdividing_sequence[i]
        end
        stage_number[n] = 8 * Int(s) - n + 3
    end
    sigma = 9 // 10

    # Initialize the constant cache
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    work = fill(zero(eltype(Q)), alg.max_order + 1)
    dt_new = fill(zero(eltype(Q)), alg.max_order + 1)
    ImplicitHairerWannerExtrapolationConstantCache(Q, n_curr, n_old, coefficients,
        stage_number, sigma, tf, uf, work,
        dt_new)
end

@cache mutable struct ImplicitHairerWannerExtrapolationCache{uType, uNoUnitsType, rateType,
    QType,
    extrapolation_coefficients,
    JType, WType, F, JCType,
    GCType, TFType, UFType} <:
                      ExtrapolationMutableCache
    # Values that are mutated
    utilde::uType
    u_temp1::uType
    u_temp2::uType
    u_temp3::Array{uType, 1}
    u_temp4::Array{uType, 1}
    tmp::uType # for get_tmp_cache()
    T::Array{uType, 1}  # Storage for the internal discretisations obtained by the explicit midpoint rule
    res::uNoUnitsType # Storage for the scaled residual of u and utilde

    fsalfirst::rateType
    k::rateType
    k_tmps::Array{rateType, 1}

    # Constant values
    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
    n_curr::Int # Storage for the current extrapolation order
    n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter
    coefficients::extrapolation_coefficients
    stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n - 1)
    sigma::Rational{Int} # Parameter for order selection

    du1::rateType
    du2::rateType
    J::JType
    W::WType
    tf::TFType
    uf::UFType
    linsolve_tmps::Array{rateType, 1}
    linsolve::Array{F, 1}
    jac_config::JCType
    grad_config::GCType
    # Values to check overflow in T1 computation
    diff1::Array{uType, 1}
    diff2::Array{uType, 1}

    #Stepsizing caches
    work::Array{QType, 1}
    dt_new::Array{QType, 1}
end

function alg_cache(alg::ImplicitHairerWannerExtrapolation, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # Initialize cache's members
    utilde = zero(u)
    u_temp1 = zero(u)
    u_temp2 = zero(u)
    u_temp3 = Array{typeof(u), 1}(undef, Threads.nthreads())
    u_temp4 = Array{typeof(u), 1}(undef, Threads.nthreads())

    for i in 1:Threads.nthreads()
        u_temp3[i] = zero(u)
        u_temp4[i] = zero(u)
    end
    tmp = zero(u)
    T = Vector{typeof(u)}(undef, alg.max_order + 1)
    for i in 1:(alg.max_order + 1)
        T[i] = zero(u)
    end
    res = uEltypeNoUnits.(zero(u))
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    k_tmps = Array{typeof(k), 1}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        k_tmps[i] = zero(rate_prototype)
    end

    cc = alg_cache(alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(false))

    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)

    if DiffEqBase.has_jac(f) && !DiffEqBase.has_Wfact(f) && f.jac_prototype !== nothing
        W_el = WOperator(f, dt, true)
        J = nothing # is J = W.J better?
    else
        J = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
        W_el = zero(J)
    end

    W = Array{typeof(W_el), 1}(undef, Threads.nthreads())
    W[1] = W_el
    for i in 2:Threads.nthreads()
        if W_el isa WOperator
            W[i] = WOperator(f, dt, true)
        else
            W[i] = zero(W_el)
        end
    end

    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)
    linsolve_tmp = zero(rate_prototype)
    linsolve_tmps = Array{typeof(linsolve_tmp), 1}(undef, Threads.nthreads())

    for i in 1:Threads.nthreads()
        linsolve_tmps[i] = zero(rate_prototype)
    end

    linprob = LinearProblem(W[1], _vec(linsolve_tmps[1]); u0 = _vec(k_tmps[1]))
    linsolve1 = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true))
    #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
    #Pr = Diagonal(_vec(weight)))

    linsolve = Array{typeof(linsolve1), 1}(undef, Threads.nthreads())
    linsolve[1] = linsolve1
    for i in 2:Threads.nthreads()
        linprob = LinearProblem(W[i], _vec(linsolve_tmps[i]); u0 = _vec(k_tmps[i]))
        linsolve[i] = init(linprob, alg.linsolve,
            alias = LinearAliasSpecifier(alias_A = true, alias_b = true))
        #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
        #Pr = Diagonal(_vec(weight)))
    end
    grad_config = build_grad_config(alg, f, tf, du1, t)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, du1, du2)

    diff1 = Array{typeof(u), 1}(undef, Threads.nthreads())
    diff2 = Array{typeof(u), 1}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        diff1[i] = zero(u)
        diff2[i] = zero(u)
    end

    # Initialize the cache
    ImplicitHairerWannerExtrapolationCache(utilde, u_temp1, u_temp2, u_temp3, u_temp4, tmp,
        T, res, fsalfirst, k, k_tmps,
        cc.Q, cc.n_curr, cc.n_old, cc.coefficients,
        cc.stage_number, cc.sigma, du1, du2, J, W, tf,
        uf, linsolve_tmps,
        linsolve, jac_config, grad_config, diff1, diff2,
        cc.work, cc.dt_new)
end

@cache mutable struct ImplicitEulerBarycentricExtrapolationConstantCache{QType,
    extrapolation_coefficients,
    TF, UF} <:
                      OrdinaryDiffEqConstantCache
    # Values that are mutated
    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
    n_curr::Int # Storage for the current extrapolation order
    n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter

    # Constant values
    coefficients::extrapolation_coefficients
    stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n - 1)
    sigma::Rational{Int} # Parameter for order selection

    tf::TF
    uf::UF

    #Stepsizing caches
    work::Array{QType, 1}
    dt_new::Array{QType, 1}
end

function alg_cache(alg::ImplicitEulerBarycentricExtrapolation, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # Initialize cache's members
    QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl

    Q = fill(zero(QType), alg.max_order + 1)
    n_curr = alg.init_order
    n_old = alg.init_order
    sequence_factor = alg.sequence_factor

    coefficients = create_extrapolation_coefficients(constvalue(uBottomEltypeNoUnits), alg)

    stage_number = Vector{Int}(undef, alg.max_order + 1)
    for n in 1:length(stage_number)
        s = zero(eltype(coefficients.subdividing_sequence))
        for i in 1:n
            s += coefficients.subdividing_sequence[i]
        end
        stage_number[n] = 2 * sequence_factor * Int(s) - n + 7
    end
    sigma = 9 // 10

    work = fill(zero(eltype(Q)), alg.max_order + 1)
    dt_new = fill(zero(eltype(Q)), alg.max_order + 1)
    # Initialize the constant cache
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    ImplicitEulerBarycentricExtrapolationConstantCache(Q, n_curr, n_old, coefficients,
        stage_number, sigma, tf, uf, work,
        dt_new)
end

@cache mutable struct ImplicitEulerBarycentricExtrapolationCache{uType, uNoUnitsType,
    rateType, QType,
    extrapolation_coefficients,
    JType, WType, F, JCType,
    GCType, TFType, UFType} <:
                      ExtrapolationMutableCache
    # Values that are mutated
    utilde::uType
    u_temp1::uType
    u_temp2::uType
    u_temp3::Array{uType, 1}
    u_temp4::Array{uType, 1}
    tmp::uType # for get_tmp_cache()
    T::Array{uType, 1}  # Storage for the internal discretisations obtained by the explicit midpoint rule
    res::uNoUnitsType # Storage for the scaled residual of u and utilde

    fsalfirst::rateType
    k::rateType
    k_tmps::Array{rateType, 1}

    # Constant values
    Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
    n_curr::Int # Storage for the current extrapolation order
    n_old::Int # Storage for the extrapolation order n_curr before perfom_step! changes the latter
    coefficients::extrapolation_coefficients
    stage_number::Vector{Int} # stage_number[n] contains information for extrapolation order (n - 1)
    sigma::Rational{Int} # Parameter for order selection

    du1::rateType
    du2::rateType
    J::JType
    W::WType
    tf::TFType
    uf::UFType
    linsolve_tmps::Array{rateType, 1}
    linsolve::Array{F, 1}
    jac_config::JCType
    grad_config::GCType
    # Values to check overflow in T1 computation
    diff1::Array{uType, 1}
    diff2::Array{uType, 1}
    #Stepsizing caches
    work::Array{QType, 1}
    dt_new::Array{QType, 1}
end

function alg_cache(alg::ImplicitEulerBarycentricExtrapolation, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # Initialize cache's members
    utilde = zero(u)
    u_temp1 = zero(u)
    u_temp2 = zero(u)
    u_temp3 = Array{typeof(u), 1}(undef, Threads.nthreads())
    u_temp4 = Array{typeof(u), 1}(undef, Threads.nthreads())

    for i in 1:Threads.nthreads()
        u_temp3[i] = zero(u)
        u_temp4[i] = zero(u)
    end
    tmp = zero(u)
    T = Vector{typeof(u)}(undef, alg.max_order + 1)
    for i in 1:(alg.max_order + 1)
        T[i] = zero(u)
    end
    res = uEltypeNoUnits.(zero(u))
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    k_tmps = Array{typeof(k), 1}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        k_tmps[i] = zero(rate_prototype)
    end

    cc = alg_cache(alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(false))

    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)

    if DiffEqBase.has_jac(f) && !DiffEqBase.has_Wfact(f) && f.jac_prototype !== nothing
        W_el = WOperator(f, dt, true)
        J = nothing # is J = W.J better?
    else
        J = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
        W_el = zero(J)
    end

    W = Array{typeof(W_el), 1}(undef, Threads.nthreads())
    W[1] = W_el
    for i in 2:Threads.nthreads()
        if W_el isa WOperator
            W[i] = WOperator(f, dt, true)
        else
            W[i] = zero(W_el)
        end
    end

    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)
    linsolve_tmp = zero(rate_prototype)
    linsolve_tmps = Array{typeof(linsolve_tmp), 1}(undef, Threads.nthreads())

    for i in 1:Threads.nthreads()
        linsolve_tmps[i] = zero(rate_prototype)
    end

    linprob = LinearProblem(W[1], _vec(linsolve_tmps[1]); u0 = _vec(k_tmps[1]))
    linsolve1 = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true))
    #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
    #Pr = Diagonal(_vec(weight)))

    linsolve = Array{typeof(linsolve1), 1}(undef, Threads.nthreads())
    linsolve[1] = linsolve1
    for i in 2:Threads.nthreads()
        linprob = LinearProblem(W[i], _vec(linsolve_tmps[i]); u0 = _vec(k_tmps[i]))
        linsolve[i] = init(linprob, alg.linsolve,
            alias = LinearAliasSpecifier(alias_A = true, alias_b = true))
        #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
        #Pr = Diagonal(_vec(weight)))
    end
    grad_config = build_grad_config(alg, f, tf, du1, t)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, du1, du2)

    diff1 = Array{typeof(u), 1}(undef, Threads.nthreads())
    diff2 = Array{typeof(u), 1}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        diff1[i] = zero(u)
        diff2[i] = zero(u)
    end

    # Initialize the cache
    ImplicitEulerBarycentricExtrapolationCache(utilde, u_temp1, u_temp2, u_temp3, u_temp4,
        tmp, T, res, fsalfirst, k, k_tmps,
        cc.Q, cc.n_curr, cc.n_old, cc.coefficients,
        cc.stage_number, cc.sigma, du1, du2, J, W,
        tf, uf, linsolve_tmps,
        linsolve, jac_config, grad_config, diff1,
        diff2, cc.work, cc.dt_new)
end
