using DelayDiffEq, DDEProblemLibrary, DiffEqDevTools
using LinearAlgebra
using Test

const prob = prob_dde_constant_2delays_long_ip
const prob_oop = prob_dde_constant_2delays_long_oop
const prob_scalar = prob_dde_constant_2delays_long_scalar

const testsol = TestSolution(
    solve(
        prob, MethodOfSteps(Vern9());
        abstol = 1 / 10^14, reltol = 1 / 10^14
    )
)

@testset "NLFunctional" begin
    alg = MethodOfSteps(Tsit5(); fpsolve = NLFunctional(; max_iter = 10))

    ## in-place problem

    sol = solve(prob, alg)

    # check statistics
    @test sol.stats.nf > 3000
    @test sol.stats.nsolve == 0
    @test sol.stats.nfpiter > 300
    @test sol.stats.nfpconvfail > 50

    # compare it with the test solution
    sol2 = appxtrue(sol, testsol)
    @test sol2.errors[:L∞] < 5.0e-4

    ## out-of-place problem

    sol_oop = solve(prob_oop, alg)

    # compare it with the in-place solution
    @test sol_oop.stats.nf == sol.stats.nf
    @test sol_oop.stats.nsolve == sol.stats.nsolve
    @test sol_oop.stats.nfpiter == sol.stats.nfpiter
    @test sol_oop.stats.nfpconvfail == sol.stats.nfpconvfail
    @test sol_oop.t ≈ sol.t atol = 1.0e-3
    @test sol_oop.u ≈ sol.u
    @test isapprox(sol.u, sol_oop.u; atol = 1.0e-7)

    ## scalar problem

    sol_scalar = solve(prob_scalar, alg)

    # compare it with the in-place solution
    @test sol_scalar.stats.nf == sol.stats.nf
    @test sol_scalar.stats.nsolve == sol.stats.nsolve
    @test sol_scalar.stats.nfpiter == sol.stats.nfpiter
    @test sol_scalar.stats.nfpconvfail == sol.stats.nfpconvfail
    @test sol_scalar.t ≈ sol.t atol = 1.0e-3
    @test sol_scalar.u ≈ sol[1, :]
end

@testset "NLAnderson" begin
    alg = MethodOfSteps(Tsit5(); fpsolve = NLAnderson(; max_iter = 10))

    ## in-place problem

    sol = solve(prob, alg)

    # check statistics
    @test sol.stats.nf < 2500
    @test sol.stats.nsolve > 0
    @test sol.stats.nfpiter < 250
    @test sol.stats.nfpconvfail < 50

    # compare it with the test solution
    sol2 = appxtrue(sol, testsol)
    @test sol2.errors[:L∞] < 5.0e-4

    ## out-of-place problem

    sol_oop = solve(prob_oop, alg)

    # compare it with the in-place solution
    @test_broken sol_oop.stats.nf == sol.stats.nf
    @test_broken sol_oop.stats.nsolve == sol.stats.nsolve
    @test_broken sol_oop.stats.nfpiter == sol.stats.nfpiter
    @test_broken sol_oop.stats.nfpconvfail == sol.stats.nfpconvfail
    @test_broken sol_oop.t ≈ sol.t
    @test_broken sol_oop.u ≈ sol.u
    @test appxtrue(sol, sol_oop).errors[:L∞] < 3.0e-6

    ## scalar problem

    sol_scalar = solve(prob_scalar, alg)

    # compare it with the in-place solution
    @test_broken sol_scalar.stats.nf == sol.stats.nf
    @test_broken sol_scalar.stats.nsolve == sol.stats.nsolve
    @test_broken sol_scalar.stats.nfpiter == sol.stats.nfpiter
    @test_broken sol_scalar.stats.nfpconvfail == sol.stats.nfpconvfail
    @test_broken sol_scalar.t ≈ sol.t
    @test_broken sol_scalar.u ≈ sol[1, :]
end
