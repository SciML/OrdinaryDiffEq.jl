using OrdinaryDiffEq, Test, Random, LinearAlgebra, SparseArrays
let N = 20
    Random.seed!(0)
    u0 = normalize(randn(N))
    dd = -2 * ones(N)
    du = ones(N - 1)
    A = diagm(-1 => du, 0 => dd, 1 => du)
    _f = (u, p, t) -> A * u - u .^ 3
    _f_ip = (du, u, p, t) -> (mul!(du, A, u); du .-= u .^ 3)
    _jac = (A, u, p, t) -> A - 3 * diagm(0 => u .^ 2)
    _jac_ip! = (J, u, p, t) -> begin
        copyto!(J, A)
        @inbounds for i in 1:N
            J[i, i] -= 3 * u[i]^2
        end
    end
    # f = ODEFunction(_f; jac=_jac)
    # f_ip = ODEFunction(_f_ip; jac=_jac_ip!, jac_prototype=zeros(N,N))
    jac_prototype = MatrixOperator(zeros(N, N); update_func! = _jac_ip!, update_func = _jac)
    f = ODEFunction(_f; jac_prototype = jac_prototype)
    f_ip = ODEFunction(_f_ip; jac_prototype = jac_prototype)
    prob = ODEProblem(f, u0, (0.0, 1.0))
    prob_ip = ODEProblem(f_ip, u0, (0.0, 1.0))

    @testset "Classical ExpRK - Low Order" begin
        dt = 0.01
        tol = 1e-3
        Algs = [LawsonEuler, NorsettEuler, ETDRK2]
        for Alg in Algs
            sol = solve(prob, Alg(krylov = true, m = 20); dt = dt, reltol = tol)
            sol_ref = solve(prob, Tsit5(); reltol = tol)
            @test_skip isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

            sol = solve(prob_ip, Alg(krylov = true, m = 20); dt = dt, reltol = tol)
            sol_ref = solve(prob_ip, Tsit5(); reltol = tol)
            @test isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

            println(Alg) # prevent Travis hanging
        end
    end

    @testset "Classical ExpRK - High Order" begin
        dt = 0.05
        tol = 1e-5
        Algs = [ETDRK3, ETDRK4, HochOst4]
        for Alg in Algs
            sol = solve(prob, Alg(krylov = true, m = 20); dt = dt, reltol = tol)
            sol_ref = solve(prob, Tsit5(); reltol = tol)
            @test_skip isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

            sol = solve(prob_ip, Alg(krylov = true, m = 20); dt = dt, reltol = tol)
            sol_ref = solve(prob_ip, Tsit5(); reltol = tol)
            @test isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

            println(Alg) # prevent Travis hanging
        end
    end

    @testset "EPIRK" begin
        dt = 0.05
        tol = 1e-5
        Algs = [Exp4, EPIRK4s3A, EPIRK4s3B, EXPRB53s3, EPIRK5P1, EPIRK5P2]
        for Alg in Algs
            sol = solve(prob, Alg(); dt = dt, reltol = tol)
            sol_ref = solve(prob, Tsit5(); reltol = tol)
            @test_broken isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

            sol = solve(prob_ip, Alg(); dt = dt, reltol = tol)
            sol_ref = solve(prob_ip, Tsit5(); reltol = tol)
            @test isapprox(sol(1.0), sol_ref(1.0); rtol = tol)
            println(Alg) # prevent Travis hanging
        end

        sol = solve(prob, EPIRK5s3(); dt = dt, reltol = tol)
        sol_ref = solve(prob, Tsit5(); reltol = tol)
        @test_broken isapprox(sol(1.0), sol_ref(1.0); rtol = tol)

        sol = solve(prob_ip, EPIRK5s3(); dt = dt, reltol = tol)
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
            for i in eachindex(sol2)
                err = maximum(abs.(sol2[i] - interp_results[i]))
                @test_skip err < tol
            end
        end

        println("Exprb32, out-of-place")
        regression_test(prob, Exprb32(m = N), 3e-4)
        println("Exprb32, inplace")
        regression_test(prob_ip, Exprb32(m = N), 3e-4)
        println("Exprb43, out-of-place")
        regression_test(prob, Exprb43(m = N), 3e-4)
        println("Exprb43, inplace")
        regression_test(prob_ip, Exprb43(m = N), 3e-4)
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
    exp_fun = ODEFunction(f;
        jac = (u, p, t) -> A,
        analytic = (u, p, t) -> exp(t * Matrix(A)) * u)
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
    exp_fun = ODEFunction(f; jac = jac, analytic = (u, p, t) -> exp(t * A) * u)
    prob = ODEProblem(exp_fun, u0, (0.0, 1.0))
    sol = solve(prob, LawsonEuler(krylov = true, m = N); dt = 0.1)
    @test sol(1.0) ≈ exp_fun.analytic(u0, nothing, 1.0)
end
