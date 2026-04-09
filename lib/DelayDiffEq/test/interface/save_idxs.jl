using DelayDiffEq
using OrdinaryDiffEqLowOrderRK
using Test

# out-of-place problem
function f_notinplace(u, h, p, t)
    return [-h(p, t - 1 / 5)[1] + u[1]; -h(p, t - 1 / 3)[2] - h(p, t - 1 / 5)[2]]
end
const prob_notinplace = DDEProblem(
    f_notinplace, ones(2), (p, t) -> zeros(2), (0.0, 100.0),
    constant_lags = [1 / 5, 1 / 3]
)

# in-place problem
function f_inplace(du, u, h, p, t)
    du[1] = -h(p, t - 1 / 5)[1] + u[1]
    return du[2] = -h(p, t - 1 / 3)[2] - h(p, t - 1 / 5)[2]
end
const prob_inplace = DDEProblem(
    f_inplace, ones(2), (p, t) -> zeros(2), (0.0, 100.0),
    constant_lags = [1 / 5, 1 / 3]
)

const alg = MethodOfSteps(BS3())

@testset "iip: $(isinplace(prob))" for prob in (prob_notinplace, prob_inplace)
    # reference solution
    dde_int = init(prob, alg)
    sol = solve!(dde_int)

    # save all components
    @testset "all components" begin
        # without keyword argument
        @testset "without keyword" begin
            # solution, solution of DDE integrator, and solution of ODE integrator
            # contain all components
            @test length(sol.u[1]) == 2
            @test length(dde_int.sol.u[1]) == 2
            @test length(dde_int.integrator.sol.u[1]) == 2

            # solution and solution of DDE integrator are equal
            @test sol.t == dde_int.sol.t
            @test sol.u == dde_int.sol.u

            ## interpolation
            @test sol(25:100, idxs = 2).u ≈ [u[1] for u in sol(25:100, idxs = [2]).u]
        end

        # with keyword argument
        @testset "with keyword" begin
            dde_int2 = init(prob, alg; save_idxs = [1, 2])
            sol2 = solve!(dde_int2)

            # solution, solution of DDE integrator, and solution of ODE integrator
            # contain all components
            @test length(sol2.u[1]) == 2
            @test length(dde_int2.sol.u[1]) == 2
            @test length(dde_int2.integrator.sol.u[1]) == 2

            # solution and solution of DDE integrator are equal
            @test sol.t == dde_int.sol.t
            @test sol.u == dde_int.sol.u

            # interpolation
            @test sol[2, :] ≈ dde_int.integrator.sol(sol.t, idxs = 2).u
        end
    end

    # save only second component
    @testset "second component" begin
        # array index
        @testset "array index" begin
            dde_int2 = init(prob, alg; save_idxs = [2])
            sol2 = solve!(dde_int2)

            # solution and solution of DDE integrator contain only second component
            @test length(sol2.u[1]) == 1
            @test length(dde_int2.sol.u[1]) == 1

            # solution of ODE integrator contains both components
            @test length(dde_int2.integrator.sol.u[1]) == 2

            # solution and solution of DDE integrator are equal
            @test sol2.t == dde_int2.sol.t
            @test sol2.u == dde_int2.sol.u

            # interpolation
            @test sol(25:100, idxs = 2).u ≈ sol2(25:100, idxs = 1).u
            @test sol(25:100, idxs = [2]).u ≈ sol2(25:100, idxs = [1]).u
        end

        # scalar index
        @testset "scalar index" begin
            dde_int2 = init(prob, alg; save_idxs = 2)
            sol2 = solve!(dde_int2)

            # solution and solution of DDE integrator is only vector of floats
            @test typeof(sol2.u) === Vector{Float64}
            @test typeof(dde_int2.sol.u) === Vector{Float64}

            # solution of ODE integrator contains both components
            @test length(dde_int2.integrator.sol.u[1]) == 2

            # solution and solution of DDE integrator are equal
            @test sol2.t == dde_int2.sol.t
            @test sol2.u == dde_int2.sol.u

            # solution equals second component of complete solution
            @test sol.t ≈ sol2.t && sol[2, :] ≈ sol2.u

            # interpolation of solution equals second component of
            # interpolation of complete solution
            @test sol(25:100, idxs = 2).u ≈ sol2(25:100, idxs = 1).u
        end
    end
end
