using StochasticDiffEq, Test

# SciML/OrdinaryDiffEq.jl#3166 / StochasticDiffEq.jl#673:
# SDE perform_step! must count drift (nf) and diffusion (nf2) evaluations.

@testset "SDE stats.nf / nf2 (#3166)" begin
    dt = 1 // 2^4
    tspan = (0.0, 1.0)
    nsteps = 16

    @testset "EM out-of-place" begin
        nf = Ref(0)
        ng = Ref(0)
        f(u, p, t) = (nf[] += 1; 1.0 * u)
        g(u, p, t) = (ng[] += 1; 1.0 * u)
        prob = SDEProblem(f, g, 0.5, tspan)
        sol = solve(prob, EM(); dt = dt)
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.naccept == nsteps
        @test sol.stats.nf == nf[] == nsteps
        @test sol.stats.nf2 == ng[] == nsteps
    end

    @testset "EM in-place" begin
        nf = Ref(0)
        ng = Ref(0)
        function f!(du, u, p, t)
            nf[] += 1
            du[1] = 1.0 * u[1]
            return nothing
        end
        function g!(du, u, p, t)
            ng[] += 1
            du[1] = 1.0 * u[1]
            return nothing
        end
        prob = SDEProblem(f!, g!, [0.5], tspan)
        sol = solve(prob, EM(); dt = dt)
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.naccept == nsteps
        @test sol.stats.nf == nf[] == nsteps
        @test sol.stats.nf2 == ng[] == nsteps
    end

    @testset "EulerHeun out-of-place" begin
        nf = Ref(0)
        ng = Ref(0)
        f(u, p, t) = (nf[] += 1; -u)
        g(u, p, t) = (ng[] += 1; 0.1 * u)
        prob = SDEProblem(f, g, 1.0, tspan)
        sol = solve(prob, EulerHeun(); dt = dt)
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.naccept == nsteps
        @test sol.stats.nf == nf[] == 2 * nsteps
        @test sol.stats.nf2 == ng[] == 2 * nsteps
    end
end
