using Test
using SimpleImplicitDiscreteSolve
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @testset "Implicit Discrete Stepping" begin
        dt = 0.01
        tspan = (0.0, 0.5)

        @testset "Linear scalar equation" begin
            function linear!(resid, u_next, u, p, t)
                # Implicit Euler for u' = -u:
                # u_{n+1} - u_n - dt*(-u_{n+1}) = 0
                resid[1] = u_next[1] - u[1] + dt * u_next[1]
                nothing
            end
            u0 = [1.0]
            nsteps = Int(round((tspan[2] - tspan[1]) / dt))
            # SimpleIDSolve initializes by solving one implicit step at t0.
            expected = u0[1] / (1 + dt)^(nsteps + 1)

            idprob = ImplicitDiscreteProblem(linear!, u0, tspan, [])
            idsol = solve(idprob, SimpleIDSolve(); dt)

            @test isapprox(idsol.u[end][1], expected, rtol = 1.0e-10)
        end

        ### free-fall
        # y, dy
        function ff(u, p, t)
            [u[2], -9.8]
        end

        function g!(resid, u_next, u, p, t)
            f = ff(u_next, p, t)
            resid[1] = u_next[1] - u[1] - 0.01 * f[1]
            resid[2] = u_next[2] - u[2] - 0.01 * f[2]
            nothing
        end
        u0 = [10.0, 0.0]
        nsteps = Int(round((tspan[2] - tspan[1]) / dt))

        idprob = ImplicitDiscreteProblem(g!, u0, tspan, [])
        idsol = solve(idprob, SimpleIDSolve(); dt)

        m = nsteps + 1
        expected_v = u0[2] - 9.8 * dt * m
        expected_y = u0[1] + m * dt * u0[2] - 9.8 * dt^2 * m * (m + 1) / 2
        @test isapprox(idsol.u[end][1], expected_y, atol = 1.0e-10)
        @test isapprox(idsol.u[end][2], expected_v, atol = 1.0e-10)
    end

    @testset "Solver initializes" begin
        function periodic!(resid, u_next, u, p, t)
            resid[1] = u_next[1] - u[1] - sin(t * π / 4)
            resid[2] = 16 - u_next[2]^2 - u_next[1]^2
        end

        tsteps = 15
        u0 = [1.0, 3.0]
        idprob = ImplicitDiscreteProblem(periodic!, u0, (0, tsteps), [])
        initsol, initfail = SciMLBase.__init(idprob, SimpleIDSolve())
        @test initsol.u[1]^2 + initsol.u[2]^2 ≈ 16

        idsol = solve(idprob, SimpleIDSolve())

        for ts in 1:tsteps
            step = idsol.u[ts]
            @test step[1]^2 + step[2]^2 ≈ 16
        end
    end
end

# Run QA tests (JET)
if TEST_GROUP != "FUNCTIONAL" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
end
