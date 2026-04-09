using DelayDiffEq, DDEProblemLibrary, Unitful
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqNonlinearSolve: NLNewton, NLFunctional, NLAnderson
using Test

using DDEProblemLibrary: remake_dde_constant_u0_tType

const probs = Dict(
    true => remake_dde_constant_u0_tType(
        prob_dde_constant_1delay_long_ip,
        [1.0u"N"],
        typeof(1.0u"s")
    ),
    false => remake_dde_constant_u0_tType(
        prob_dde_constant_1delay_long_scalar,
        1.0u"N",
        typeof(1.0u"s")
    )
)

# we test the current handling of units for regressions
# however, it is broken upstream: https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/issues/828

@testset "iip: $inplace" for inplace in (true, false)
    prob = probs[inplace]

    alg = MethodOfSteps(
        Tsit5(); constrained = false,
        fpsolve = NLFunctional(; max_iter = 100)
    )

    # default
    sol1 = solve(prob, alg)

    # without units
    if inplace
        @test_throws Unitful.DimensionError solve(
            prob, alg;
            abstol = 1.0e-6, reltol = 1.0e-3
        )
    else
        sol2 = solve(prob, alg; abstol = 1.0e-6, reltol = 1.0e-3)

        @test sol1.t == sol2.t
        @test sol1.u == sol2.u
    end

    # with correct units
    sol3 = solve(prob, alg; abstol = 1.0e-6u"N", reltol = 1.0e-3u"N")

    @test sol1.t == sol3.t
    @test sol1.u == sol3.u

    # with correct units as vectors
    if inplace
        sol4 = solve(prob, alg; abstol = [1.0e-6u"N"], reltol = [1.0e-3u"N"])

        @test sol1.t == sol4.t
        @test sol1.u == sol4.u
    end

    # with incorrect units for absolute tolerance
    @test_throws Unitful.DimensionError solve(
        prob, alg;
        abstol = 1.0e-6u"s", reltol = 1.0e-3u"N"
    )

    # with incorrect units for relative tolerance
    @test_throws Unitful.DimensionError solve(
        prob, alg;
        abstol = 1.0e-6u"N", reltol = 1.0e-3u"s"
    )
end
