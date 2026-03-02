using DelayDiffEq
using Test

@testset "in-place" begin
    # define functions (Hutchinson's equation)
    function f(du, u, h, p, t)
        du[1] = u[1] * (1 - h(p, t - 1)[1])
        nothing
    end

    njacs = Ref(0)
    function jac(J, u, h, p, t)
        njacs[] += 1
        J[1, 1] = 1 - h(p, t - 1)[1]
        nothing
    end

    nWfact_ts = Ref(0)
    function Wfact_t(W, u, h, p, dtgamma, t)
        nWfact_ts[] += 1
        W[1, 1] = 1 - h(p, t - 1)[1] - inv(dtgamma)
        nothing
    end

    h(p, t) = [0.0]

    # define problems
    prob = DDEProblem(DDEFunction{true}(f), [1.0], h, (0.0, 40.0); constant_lags = [1])
    prob_jac = remake(prob; f = DDEFunction{true}(f; jac = jac))
    prob_Wfact_t = remake(prob; f = DDEFunction{true}(f; Wfact_t = Wfact_t))

    # compute solutions
    for alg in (Rodas5P(), TRBDF2())
        sol = solve(prob, MethodOfSteps(alg))

        ## Jacobian
        njacs[] = 0
        sol_jac = solve(prob_jac, MethodOfSteps(alg))

        # check number of function evaluations
        @test !iszero(njacs[])
        @test njacs[] == sol_jac.stats.njacs
        @test njacs[] <= sol_jac.stats.nw

        # check resulting solution
        @test sol.t ≈ sol_jac.t
        @test sol.u ≈ sol_jac.u

        ## Wfact_t
        nWfact_ts[] = 0
        sol_Wfact_t = solve(prob_Wfact_t, MethodOfSteps(alg))

        # check number of function evaluations
        @test !iszero(nWfact_ts[])
        @test_broken nWfact_ts[] == njacs[]
        @test iszero(sol_Wfact_t.stats.njacs)
        @test_broken nWfact_ts[] == sol_Wfact_t.stats.nw

        # check resulting solution
        if alg isa TRBDF2
            @test_broken sol.t ≈ sol_Wfact_t.t
            @test_broken sol.u ≈ sol_Wfact_t.u
        else
            @test sol.t ≈ sol_Wfact_t.t
            @test sol.u ≈ sol_Wfact_t.u
        end
    end
end

@testset "out-of-place" begin
    # define functions (Hutchinson's equation)
    f(u, h, p, t) = u[1] .* (1 .- h(p, t - 1))

    njacs = Ref(0)
    function jac(u, h, p, t)
        njacs[] += 1
        reshape(1 .- h(p, t - 1), 1, 1)
    end

    nWfact_ts = Ref(0)
    function Wfact_t(u, h, p, dtgamma, t)
        nWfact_ts[] += 1
        reshape((1 - inv(dtgamma)) .- h(p, t - 1), 1, 1)
    end

    h(p, t) = [0.0]

    # define problems
    prob = DDEProblem(DDEFunction{false}(f), [1.0], h, (0.0, 40.0); constant_lags = [1])
    prob_jac = remake(prob; f = DDEFunction{false}(f; jac = jac))
    prob_Wfact_t = remake(prob; f = DDEFunction{false}(f; Wfact_t = Wfact_t))

    # compute solutions
    # Only test Rodas5P: OOP TRBDF2 hangs with OrdinaryDiffEqCore >= 3.10
    for alg in (Rodas5P(),)
        sol = solve(prob, MethodOfSteps(alg))

        ## Jacobian
        njacs[] = 0
        sol_jac = solve(prob_jac, MethodOfSteps(alg))

        # check number of function evaluations
        @test !iszero(njacs[])
        @test_broken njacs[] == sol_jac.stats.njacs
        @test_broken njacs[] == sol_jac.stats.nw

        # check resulting solution
        @test sol.t ≈ sol_jac.t
        @test sol.u ≈ sol_jac.u

        ## Wfact_t
        nWfact_ts[] = 0
        sol_Wfact_t = solve(prob_Wfact_t, MethodOfSteps(alg))

        # check number of function evaluations
        @test_broken !iszero(nWfact_ts[])
        @test_broken nWfact_ts[] == njacs[]
        @test_broken iszero(sol_Wfact_ts.stats.njacs)
        @test_broken nWfact_ts[] == sol_Wfact_t.stats.nw

        # check resulting solution
        @test sol.t ≈ sol_Wfact_t.t
        @test sol.u ≈ sol_Wfact_t.u
    end
end
