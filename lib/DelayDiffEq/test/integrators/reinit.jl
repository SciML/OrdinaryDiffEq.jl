using DelayDiffEq, DDEProblemLibrary, RecursiveArrayTools
using OrdinaryDiffEqLowOrderRK
using Test

const prob_ip = prob_dde_constant_1delay_ip
const prob_scalar = prob_dde_constant_1delay_scalar

const alg = MethodOfSteps(BS3(); constrained = false)

@testset "iip: $inplace" for inplace in (true, false)
    prob = inplace ? prob_ip : prob_scalar

    @testset "integrator" begin
        integrator = init(prob, alg, dt = 0.01)
        solve!(integrator)

        u = recursivecopy(integrator.sol.u)
        t = copy(integrator.sol.t)

        reinit!(integrator)
        integrator.dt = 0.01
        solve!(integrator)

        @test u == integrator.sol.u
        @test t == integrator.sol.t
    end

    @testset "solution" begin
        integrator = init(prob, alg, dt = 0.01, tstops = [0.5], saveat = [0.33])
        sol = solve!(integrator)

        u = recursivecopy(sol.u)
        t = copy(sol.t)

        reinit!(integrator)
        integrator.dt = 0.01
        sol = solve!(integrator)

        @test u == sol.u
        @test t == sol.t
    end

    # https://github.com/SciML/OrdinaryDiffEq.jl/issues/1394
    @testset "saveiter_dense" begin
        integrator = init(prob, alg)
        solve!(integrator)

        # ensure that the counters were updated
        for i in (integrator, integrator.integrator)
            @test i.saveiter > 1
            @test i.saveiter_dense > 1
        end

        reinit!(integrator)

        # check that the counters are reset
        for i in (integrator, integrator.integrator)
            @test i.saveiter == 1
            @test i.saveiter_dense == 1
        end
    end
end
