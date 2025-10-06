using Test
using SimpleImplicitDiscreteSolve
using OrdinaryDiffEqSDIRK

# Test implicit Euler using ImplicitDiscreteProblem
@testset "Implicit Euler" begin
    function lotkavolterra(u, p, t)
        [1.5 * u[1] - u[1] * u[2], -3.0 * u[2] + u[1] * u[2]]
    end

    function f!(resid, u_next, u, p, t)
        lv = lotkavolterra(u_next, p, t)
        resid[1] = u_next[1] - u[1] - 0.01 * lv[1]
        resid[2] = u_next[2] - u[2] - 0.01 * lv[2]
        nothing
    end
    u0 = [1.0, 1.0]
    tspan = (0.0, 0.5)

    idprob = ImplicitDiscreteProblem(f!, u0, tspan, [])
    idsol = solve(idprob, SimpleIDSolve(), dt = 0.01)

    oprob = ODEProblem(lotkavolterra, u0, tspan)
    osol = solve(oprob, ImplicitEuler())

    @test isapprox(idsol[end - 1], osol[end], atol = 0.1)

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
    tspan = (0, 0.5)

    idprob = ImplicitDiscreteProblem(g!, u0, tspan, [])
    idsol = solve(idprob, SimpleIDSolve(); dt = 0.01)

    oprob = ODEProblem(ff, u0, tspan)
    osol = solve(oprob, ImplicitEuler())

    @test isapprox(idsol[end - 1], osol[end], atol = 0.1)
end

@testset "Solver initializes" begin
    function periodic!(resid, u_next, u, p, t)
        resid[1] = u_next[1] - u[1] - sin(t * π / 4)
        resid[2] = 16 - u_next[2]^2 - u_next[1]^2
    end

    tsteps = 15
    u0 = [1.0, 3.0]
    idprob = ImplicitDiscreteProblem(periodic!, u0, (0, tsteps), [])
    initsol, initfail = DiffEqBase.__init(idprob, SimpleIDSolve())
    @test initsol.u[1]^2 + initsol.u[2]^2 ≈ 16

    idsol = solve(idprob, SimpleIDSolve())

    for ts in 1:tsteps
        step = idsol.u[ts]
        @test step[1]^2 + step[2]^2 ≈ 16
    end
end
