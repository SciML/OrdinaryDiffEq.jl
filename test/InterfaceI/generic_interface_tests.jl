import CommonSolve
import SciMLBase
using OrdinaryDiffEqTsit5: Tsit5
using Test

@testset "generic SciML solve interface" begin
    rhs(u, p, t) = p * u
    prob = SciMLBase.ODEProblem(rhs, 1.0, (0.0, 1.0), -1.0)

    @testset "out-of-place solve and remake" begin
        sol = CommonSolve.solve(prob, Tsit5(); saveat = [0.0, 0.5, 1.0])
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end] ≈ exp(-1) rtol = 1.0e-6
        @test sol.t == [0.0, 0.5, 1.0]

        remade = SciMLBase.remake(prob; p = -2.0, tspan = (0.0, 0.5))
        remade_sol = CommonSolve.solve(remade, Tsit5())
        @test SciMLBase.successful_retcode(remade_sol)
        @test remade_sol.u[end] ≈ exp(-1) rtol = 1.0e-6
    end

    @testset "in-place solve" begin
        function rhs!(du, u, p, t)
            du[1] = p * u[1]
            return nothing
        end

        inplace_prob = SciMLBase.ODEProblem(rhs!, [1.0], (0.0, 1.0), -1.0)
        sol = CommonSolve.solve(inplace_prob, Tsit5())
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end][1] ≈ exp(-1) rtol = 1.0e-6
    end

    @testset "generic integrator lifecycle and time stops" begin
        integrator = CommonSolve.init(
            prob,
            Tsit5();
            adaptive = false,
            dt = 0.1,
            save_everystep = true,
        )
        SciMLBase.add_tstop!(integrator, 0.25)
        SciMLBase.step!(integrator)
        @test integrator.t > 0
        CommonSolve.solve!(integrator)
        @test SciMLBase.successful_retcode(integrator.sol)
        @test 0.25 ∈ integrator.sol.t
    end
end
