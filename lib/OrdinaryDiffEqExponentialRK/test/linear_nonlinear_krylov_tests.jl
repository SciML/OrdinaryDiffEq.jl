using OrdinaryDiffEqExponentialRK, Test, Random, LinearAlgebra, SparseArrays
using OrdinaryDiffEqTsit5
using SciMLOperators: MatrixOperator, isconstant
using SciMLBase: successful_retcode
using OrdinaryDiffEqExponentialRK: _cached_ishermitian, _arnoldi_kwargs

let N = 20
    Random.seed!(0)
    u0 = normalize(randn(N))
    dd = -2 * ones(N)
    du = ones(N - 1)
    A = diagm(-1 => du, 0 => dd, 1 => du)
    _f = (u, p, t) -> A * u - u .^ 3
    _f_ip = (du, u, p, t) -> (mul!(du, A, u); du .-= u .^ 3)
    _jac = (J, u, p, t) -> A - 3 * diagm(0 => u .^ 2)
    _jac_ip! = (J, u, p, t) -> begin
        copyto!(J, A)
        @inbounds for i in 1:N
            J[i, i] -= 3 * u[i]^2
        end
    end
    # f = ODEFunction(_f; jac=_jac)
    # f_ip = ODEFunction(_f_ip; jac=_jac_ip!, jac_prototype=zeros(N,N))
    jac_prototype = MatrixOperator(zeros(N, N); update_func! = _jac_ip!, update_func = _jac)
    f = ODEFunction(_f; jac_prototype)
    f_ip = ODEFunction(_f_ip; jac_prototype)
    prob = ODEProblem(f, u0, (0.0, 1.0))
    prob_ip = ODEProblem(f_ip, u0, (0.0, 1.0))

    @testset "Classical ExpRK - Low Order" begin
        dt = 0.01
        tol = 1.0e-3
        Algs = [LawsonEuler, NorsettEuler, ETDRK2]
        for Alg in Algs
            sol = solve(prob, Alg(krylov = true, m = 20); dt, reltol = tol)
            sol_ref = solve(prob, Tsit5(); reltol = tol)
            @test isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

            sol = solve(prob_ip, Alg(krylov = true, m = 20); dt, reltol = tol)
            sol_ref = solve(prob_ip, Tsit5(); reltol = tol)
            @test isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

            println(Alg) # prevent Travis hanging
        end
    end

    @testset "Classical ExpRK - High Order" begin
        dt = 0.05
        tol = 1.0e-5
        Algs = [ETDRK3, ETDRK4, HochOst4, Friedli]
        for Alg in Algs
            sol = solve(prob, Alg(krylov = true, m = 20); dt, reltol = tol)
            sol_ref = solve(prob, Tsit5(); reltol = tol)
            @test isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

            sol = solve(prob_ip, Alg(krylov = true, m = 20); dt, reltol = tol)
            sol_ref = solve(prob_ip, Tsit5(); reltol = tol)
            @test isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

            println(Alg) # prevent Travis hanging
        end
    end

    @testset "EPIRK" begin
        dt = 0.05
        tol = 1.0e-5
        Algs = [Exp4, EPIRK4s3A, EPIRK4s3B, EXPRB53s3, EPIRK5P1, EPIRK5P2]
        for Alg in Algs
            sol = solve(prob, Alg(); dt, reltol = tol)
            sol_ref = solve(prob, Tsit5(); reltol = tol)
            @test isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

            sol = solve(prob_ip, Alg(); dt, reltol = tol)
            sol_ref = solve(prob_ip, Tsit5(); reltol = tol)
            @test isapprox(sol(1.0), sol_ref(1.0); rtol = tol)
            println(Alg) # prevent Travis hanging
        end

        sol = solve(prob, EPIRK5s3(); dt, reltol = tol)
        sol_ref = solve(prob, Tsit5(); reltol = tol)
        @test_broken isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

        sol = solve(prob_ip, EPIRK5s3(); dt, reltol = tol)
        sol_ref = solve(prob_ip, Tsit5(); reltol = tol)
        @test_broken isapprox(sol(1.0), sol_ref(1.0); rtol = tol)
        println(EPIRK5s3) # prevent Travis hanging
    end

    @testset "Adaptive exponential Rosenbrock" begin
        # Regression tests adapted from ode_dense_tests.jl
        interp_points = 0.0:(1 / 16):1.0
        interp_results = [zeros(N) for _ in 1:length(interp_points)]
        function regression_test(prob, alg, tol)
            sol1 = solve(prob, alg, dt = 1 / 4, dense = true, adaptive = true)
            sol1(interp_results, interp_points)
            sol2 = solve(prob, alg, dt = 1 / 16, dense = true, adaptive = false)
            for i in eachindex(interp_results)
                err = maximum(abs.(sol2.u[i] - interp_results[i]))
                @test err < tol
            end
        end

        println("Exprb32, out-of-place")
        regression_test(prob, Exprb32(m = N), 3.0e-4)
        println("Exprb32, inplace")
        regression_test(prob_ip, Exprb32(m = N), 3.0e-4)
        println("Exprb43, out-of-place")
        regression_test(prob, Exprb43(m = N), 3.0e-4)
        println("Exprb43, inplace")
        regression_test(prob_ip, Exprb43(m = N), 3.0e-4)
    end
end

@testset "ExpRK with custom jacobian" begin
    N = 10
    # Sparse Jacobian
    Random.seed!(0)
    u0 = normalize(randn(N))
    dd = -2 * ones(N)
    du = ones(N - 1)
    A = spdiagm(-1 => du, 0 => dd, 1 => du)
    f = (u, p, t) -> A * u
    exp_fun = ODEFunction(
        f;
        jac = (u, p, t) -> A,
        analytic = (u, p, t) -> exp(t * Matrix(A)) * u
    )
    prob = ODEProblem(exp_fun, u0, (0.0, 1.0))
    sol = solve(prob, LawsonEuler(krylov = true, m = N); dt = 0.1)
    @test sol(1.0) ≈ exp_fun.analytic(u0, nothing, 1.0)
end

@testset "ExpRK with default jacobian" begin
    N = 10
    Random.seed!(0)
    u0 = normalize(randn(N))
    dd = -2 * ones(N)
    du = ones(N - 1)
    A = diagm(-1 => du, 0 => dd, 1 => du)
    f = (du, u, p, t) -> mul!(du, A, u)
    jac = (J, u, p, t) -> (J .= A; nothing)
    exp_fun = ODEFunction(f; jac, analytic = (u, p, t) -> exp(t * A) * u)
    prob = ODEProblem(exp_fun, u0, (0.0, 1.0))
    sol = solve(prob, LawsonEuler(krylov = true, m = N); dt = 0.1)
    @test sol(1.0) ≈ exp_fun.analytic(u0, nothing, 1.0)
end

# Counts how often the solver asks the linear part for its symmetry. `MatrixOperator` forwards
# `ishermitian` to the wrapped array, so this sees exactly the calls `arnoldi!` makes when it is
# left to derive the flag itself.
mutable struct SymmetryCounter{T} <: AbstractMatrix{T}
    A::SparseMatrixCSC{T, Int}
    count::Int
end
SymmetryCounter(A::SparseMatrixCSC{T, Int}) where {T} = SymmetryCounter{T}(A, 0)
Base.size(C::SymmetryCounter) = size(C.A)
Base.getindex(C::SymmetryCounter, i::Int, j::Int) = C.A[i, j]
LinearAlgebra.mul!(y::AbstractVector, C::SymmetryCounter, x::AbstractVector) = mul!(y, C.A, x)
function LinearAlgebra.mul!(y::AbstractVector, C::SymmetryCounter, x::AbstractVector, a, b)
    return mul!(y, C.A, x, a, b)
end
function LinearAlgebra.ishermitian(C::SymmetryCounter)
    C.count += 1
    return ishermitian(C.A)
end

@testset "Cached ishermitian flag" begin
    # `arnoldi!` takes `ishermitian` as a default keyword argument, so it re-derives the
    # property on every call -- five times per ETDRK4 step, and for a symmetric sparse
    # operator that is a full O(nnz) scan. A *constant* linear part holds still for the whole
    # solve, so `alg_cache_expRK` evaluates it once and `_arnoldi_kwargs` threads it through.
    # These problems are split, unlike the ODEProblem fixtures above.
    Random.seed!(0)
    N = 32
    g!(du, u, p, t) = (@. du = u - u^3)
    split_prob(A) = SplitODEProblem(
        MatrixOperator(sparse(A)), g!,
        normalize(randn(N)), (0.0, 0.1)
    )

    sym = split_prob(
        begin                       # periodic Laplacian => symmetric
            A = diagm(-1 => ones(N - 1), 0 => -2ones(N), 1 => ones(N - 1)) .* 50.0
            A[1, N] = A[N, 1] = 50.0
            A
        end
    )
    nonsym = split_prob(                         # upwind-biased => not symmetric
        diagm(-1 => 3ones(N - 1), 0 => -4ones(N), 1 => ones(N - 1)) .* 50.0
    )

    # the three distinct outcomes
    @test _cached_ishermitian(sym.f) === true
    @test _cached_ishermitian(nonsym.f) === false
    @test _cached_ishermitian(
        ODEProblem(
            (du, u, p, t) -> (@. du = -u), [1.0],
            (0.0, 1.0)
        ).f
    ) === nothing

    # An operator may be symmetric at t=0 and not afterwards. A Jacobian cannot reach here
    # (the SplitFunction guard above sends it down the uncached path), but a SplitFunction
    # whose linear part carries an `update_func` can have its entries replaced by
    # `update_coefficients!`. Caching `true` there would send `arnoldi!` to `lanczos!`,
    # whose three-term recurrence assumes symmetry -- silently wrong, not merely slow. So
    # anything that is not `isconstant` must decline to cache and fall back per call.
    varying = SplitODEProblem(
        MatrixOperator(
            sparse(sym.f.f1.f.A);
            update_func = (A, u, p, t) -> (B = copy(A); B[1, 2] += t; B)
        ),
        g!, normalize(randn(N)), (0.0, 0.1)
    )

    @test ishermitian(varying.f.f1.f) === true    # symmetric at construction ...
    @test isconstant(varying.f.f1.f) === false    # ... but free to stop being so
    @test _cached_ishermitian(varying.f) === nothing

    # having declined, the per-call fallback must still produce the right answer
    @test _arnoldi_kwargs(
        ETDRK4(krylov = true, m = 15), varying.f.f1.f,
        (; opts = (; internalopnorm = opnorm)), nothing
    ).ishermitian ===
        ishermitian(varying.f.f1.f)

    # and the cache must actually hold `nothing`, not a stale snapshot
    @test init(varying, ETDRK4(krylov = true, m = 15); dt = 1.0e-3).cache.KsCache[4] ===
        nothing

    # The operator really does stop being symmetric part-way through the solve -- evaluating
    # the split right-hand side runs `update_coefficients!` on it -- so a snapshot taken at
    # t=0 would be wrong, not merely stale. What reaches `arnoldi!` has to track the operator
    # rather than the shape of the problem.
    let alg = ETDRK4(krylov = true, m = 15), integ = init(varying, alg; dt = 1.0e-3)
        step!(integ)
        step!(integ)
        A = integ.f.f1.f
        @test ishermitian(A) === false
        @test _arnoldi_kwargs(alg, A, integ, integ.cache.KsCache[4]).ishermitian === false
    end

    # when nothing was cached, the fallback must derive the same value, with the same
    # NamedTuple shape so the call sites stay type-stable
    let A = sym.f.f1.f, integ = (; opts = (; internalopnorm = opnorm)),
            alg = ETDRK4(krylov = true, m = 15)

        @test _arnoldi_kwargs(alg, A, integ, true).ishermitian ===
            _arnoldi_kwargs(alg, A, integ, nothing).ishermitian === true
        @test keys(_arnoldi_kwargs(alg, A, integ, true)) ===
            keys(_arnoldi_kwargs(alg, A, integ, nothing))
    end

    # The flag must reach every cache. Checked per algorithm rather than on one
    # representative, since each `perform_step!` builds its own `arnoldi!` keywords.
    # NOTE this asserts only that the flag is STORED: `alg_cache_expRK` populates
    # `KsCache[4]` identically for all of them, so it cannot catch a `perform_step!` that
    # receives the flag and then ignores it. All four in-place caches route through
    # `_arnoldi_kwargs`; a new scheme that hand-rolls its keywords would pass this and
    # still silently re-derive `ishermitian` on every build.
    for prob in (sym, nonsym), Alg in (ETDRK2, ETDRK3, ETDRK4, HochOst4)
        cache = init(prob, Alg(krylov = true, m = 15); dt = 1.0e-3).cache
        @test cache.KsCache[4] == ishermitian(prob.f.f1.f)
    end

    # That loop asserts only that the flag is stored. This one counts what the operator is
    # actually asked for: the symmetry check happens while the cache is built and never again,
    # however many steps run. A `perform_step!` that drops the flag re-derives it per build.
    for Alg in (NorsettEuler, ETDRK2, ETDRK3, ETDRK4, HochOst4)
        counter = SymmetryCounter(sparse(sym.f.f1.f.A))
        integ = init(
            SplitODEProblem(MatrixOperator(counter), g!, normalize(randn(N)), (0.0, 0.1)),
            Alg(krylov = true, m = 15); dt = 1.0e-3
        )
        step!(integ)
        built = counter.count
        for _ in 1:4
            step!(integ)
        end
        @test built > 0
        @test counter.count == built
    end

    # smoke: the split path still integrates (the fixtures above are all non-split)
    for prob in (sym, nonsym)
        @test successful_retcode(
            solve(
                prob, ETDRK4(krylov = true, m = 15);
                dt = 1.0e-3, save_everystep = false
            )
        )
    end
end
