using DelayDiffEq, DDEProblemLibrary
using Test

# Check that numerical solutions approximate analytical solutions,
# independent of problem structure

@testset "standard history" begin
    # standard algorithm
    alg = MethodOfSteps(BS3(); constrained = false)

    alg1 = MethodOfSteps(
        Tsit5(); constrained = false,
        fpsolve = NLFunctional(; max_iter = 100)
    )
    alg2 = MethodOfSteps(
        DP8(); constrained = false,
        fpsolve = NLFunctional(; max_iter = 10)
    )
    alg3 = MethodOfSteps(Tsit5(); constrained = true)
    alg4 = MethodOfSteps(
        DP5(); constrained = false,
        fpsolve = NLFunctional(; max_iter = 100)
    )

    ## Single constant delay
    @testset "single constant delay" begin
        @testset "short time span" begin
            ### Scalar function
            sol_scalar = solve(prob_dde_constant_1delay_scalar, alg)

            @test sol_scalar.errors[:l∞] < 3.0e-5
            @test sol_scalar.errors[:final] < 2.1e-5
            @test sol_scalar.errors[:l2] < 1.3e-5

            ### Out-of-place function
            sol_oop = solve(prob_dde_constant_1delay_oop, alg)

            @test sol_scalar.t ≈ sol_oop.t && sol_scalar.u ≈ sol_oop[1, :]

            ### In-place function
            sol_ip = solve(prob_dde_constant_1delay_ip, alg)

            @test sol_scalar.t ≈ sol_ip.t && sol_scalar.u ≈ sol_ip[1, :]
        end

        @testset "long time span" begin
            prob = prob_dde_constant_1delay_long_scalar

            sol1 = solve(prob, alg1; abstol = 1.0e-12, reltol = 1.0e-12)
            sol2 = solve(prob, alg2; abstol = 1.0e-8, reltol = 1.0e-10)
            sol3 = solve(prob, alg3; abstol = 1.0e-8, reltol = 1.0e-10)
            sol4 = solve(prob, alg4; abstol = 1.0e-12, reltol = 1.0e-12)

            # relaxed tests to prevent floating point issues
            @test abs(sol1.u[end] - sol2.u[end]) < 2.5e-8
            @test abs(sol1.u[end] - sol3.u[end]) < 3.7e-8
            @test abs(sol1.u[end] - sol4.u[end]) < 9.0e-11 # 9.0e-13
        end
    end

    ## Two constant delays
    @testset "two constant delays" begin
        @testset "short time span" begin
            ### Scalar function
            sol_scalar = solve(prob_dde_constant_2delays_scalar, alg)

            @test sol_scalar.errors[:l∞] < 2.4e-6
            @test sol_scalar.errors[:final] < 2.1e-6
            @test sol_scalar.errors[:l2] < 1.2e-6

            ### Out-of-place function
            sol_oop = solve(prob_dde_constant_2delays_oop, alg)

            @test sol_scalar.t ≈ sol_oop.t && sol_scalar.u ≈ sol_oop[1, :]

            ### In-place function
            sol_ip = solve(prob_dde_constant_2delays_ip, alg)

            @test sol_scalar.t ≈ sol_ip.t && sol_scalar.u ≈ sol_ip[1, :]
        end

        @testset "long time span" begin
            prob = prob_dde_constant_2delays_long_scalar

            sol1 = solve(prob, alg1; abstol = 1.0e-12, reltol = 1.0e-12)
            sol2 = solve(prob, alg2; abstol = 1.0e-8, reltol = 1.0e-10)
            sol3 = solve(prob, alg3; abstol = 1.0e-8, reltol = 1.0e-10)
            sol4 = solve(prob, alg4; abstol = 1.0e-12, reltol = 1.0e-12)

            # relaxed tests to prevent floating point issues
            @test abs(sol1.u[end] - sol3.u[end]) < 1.2e-13 # 1.2e-15
            @test abs(sol1.u[end] - sol4.u[end]) < 3.1e-13 # 3.1e-15
        end
    end
end

## Non-standard history functions
@testset "non-standard history" begin
    alg = MethodOfSteps(
        Tsit5(); constrained = false,
        fpsolve = NLFunctional(; max_iter = 100)
    )

    @testset "idxs" begin
        function f(du, u, h, p, t)
            du[1] = -h(p, t - 0.2; idxs = 1) + u[1]
        end
        h(p, t; idxs = nothing) = idxs isa Number ? 0.0 : [0.0]

        prob = DDEProblem(f, [1.0], h, (0.0, 100.0); constant_lags = [0.2])
        solve(prob, alg; abstol = 1.0e-12, reltol = 1.0e-12)
    end

    @testset "in-place" begin
        function f(du, u, h, p, t)
            h(du, p, t - 0.2)
            du[1] = -du[1] + u[1]
        end
        h(val, p, t) = (val .= 0.0)

        prob = DDEProblem(f, [1.0], h, (0.0, 100.0); constant_lags = [0.2])
        solve(prob, alg; abstol = 1.0e-12, reltol = 1.0e-12)
    end
end
