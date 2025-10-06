#runtests
using Test
using ImplicitDiscreteSolve
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using SciMLBase
using JET

# Test implicit Euler using ImplicitDiscreteProblem
@testset "Implicit Euler" begin
    function lotkavolterra(u, p, t)
        [1.5*u[1] - u[1]*u[2], -3.0*u[2] + u[1]*u[2]]
    end

    function f!(resid, u_next, u, p, t)
        lv = lotkavolterra(u_next, p, t)
        resid[1] = u_next[1] - u[1] - 0.01*lv[1]
        resid[2] = u_next[2] - u[2] - 0.01*lv[2]
        nothing
    end
    u0 = [1., 1.]
    tspan = (0., 0.5)

    idprob = ImplicitDiscreteProblem(f!, u0, tspan, []; dt = 0.01)
    idsol = solve(idprob, IDSolve())

    oprob = ODEProblem(lotkavolterra, u0, tspan)
    osol = solve(oprob, ImplicitEuler())

    @test isapprox(idsol[end], osol[end], atol = 0.1)

    ### free-fall
    # y, dy
    function ff(u, p, t) 
        [u[2], -9.8]
    end

    function g!(resid, u_next, u, p, t) 
        f = ff(u_next, p, t) 
        resid[1] = u_next[1] - u[1] - 0.01*f[1]
        resid[2] = u_next[2] - u[2] - 0.01*f[2]
        nothing
    end
    u0 = [10., 0.]
    tspan = (0, 0.2)

    idprob = ImplicitDiscreteProblem(g!, u0, tspan, []; dt = 0.01)
    idsol = solve(idprob, IDSolve())

    oprob = ODEProblem(ff, u0, tspan)
    osol = solve(oprob, ImplicitEuler())

    @test isapprox(idsol[end], osol[end], atol = 0.1)
end

@testset "Solver initializes" begin
    function periodic!(resid, u_next, u, p, t) 
        resid[1] = u_next[1] - u[1] - sin(t*π/4)
        resid[2] = 16 - u_next[2]^2 - u_next[1]^2 
    end

    tsteps = 15
    u0 = [1., 3.]
    idprob = ImplicitDiscreteProblem(periodic!, u0, (0, tsteps), [])
    integ = init(idprob, IDSolve())
    @test integ.u[1]^2 + integ.u[2]^2 ≈ 16

    for ts in 1:tsteps
        step!(integ)
        @test integ.u[1]^2 + integ.u[2]^2 ≈ 16
    end
end

@testset "Handle nothing in u0" begin
    function empty(u_next, u, p, t) 
        nothing
    end

    tsteps = 5
    u0 = nothing
    idprob = ImplicitDiscreteProblem(empty, u0, (0, tsteps), [])
    @test_nowarn integ = init(idprob, IDSolve())
end

@testset "Create NonlinearLeastSquaresProblem" begin
    function over(u_next, u, p, t) 
        [u_next[1] - 1, u_next[2] - 1, u_next[1] - u_next[2]]
    end

    tsteps = 5
    u0 = [1., 1.]
    idprob = ImplicitDiscreteProblem(ImplicitDiscreteFunction(over, resid_prototype = zeros(3)), u0, (0, tsteps), [])
    integ = init(idprob, IDSolve())
    @test integ.cache.prob isa NonlinearLeastSquaresProblem

    function under(u_next, u, p, t) 
        [u_next[1] - u_next[2] - 1] 
    end
    idprob = ImplicitDiscreteProblem(ImplicitDiscreteFunction(under; resid_prototype = zeros(1)), u0, (0, tsteps), [])
    integ = init(idprob, IDSolve())
    @test integ.cache.prob isa NonlinearLeastSquaresProblem

    function full(u_next, u, p, t) 
        [u_next[1]^2 - 3, u_next[2] - u[1]] 
    end
    idprob = ImplicitDiscreteProblem(ImplicitDiscreteFunction(full; resid_prototype = zeros(2)), u0, (0, tsteps), [])
    integ = init(idprob, IDSolve())
    @test integ.cache.prob isa NonlinearProblem
end

@testset "InitialFailure thrown" begin
    function bad(u_next, u, p, t) 
        [u_next[1] - u_next[2], u_next[1] - 3, u_next[2] - 4]
    end

    u0 = [3., 4.]
    idprob = ImplicitDiscreteProblem(bad, u0, (0, 0), [])
    integ = init(idprob, IDSolve())
    @test check_error(integ) == ReturnCode.InitialFailure
    sol = solve(idprob, IDSolve())
    @test length(sol.u) == 1
end


@testset "JET Tests" begin
    test_package(test_package(
        ImplicitDiscreteSolve, target_defined_modules = true, mode = :typo))
end

include("qa_tests.jl")