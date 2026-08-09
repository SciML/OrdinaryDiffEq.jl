# This file must not load RecursiveFactorization: the regression it guards
# (issue #4154) is only visible when RecursiveFactorization is absent.
using OrdinaryDiffEqExtrapolation, Polyester, Test
using CommonSolve: solve
using SciMLBase: SciMLBase, ODEProblem
using LinearSolve: GenericLUFactorization
using OrdinaryDiffEqCore: PolyesterThreads

function vanderpol!(du, u, p, t)
    du[1] = u[2]
    du[2] = p * (1 - u[1]^2) * u[2] - u[1]
    return nothing
end
prob = ODEProblem(vanderpol!, [2.0, 0.0], (0.0, 1.0), 100.0)

algs = [
    ImplicitEulerExtrapolation, ImplicitDeuflhardExtrapolation,
    ImplicitHairerWannerExtrapolation, ImplicitEulerBarycentricExtrapolation,
]

@testset "Threaded implicit extrapolation defaults to a single-threaded factorization" begin
    @testset "$Alg" for Alg in algs
        serial = Alg(threading = false)
        @test serial.linsolve === nothing
        reference = solve(prob, serial; abstol = 1.0e-9, reltol = 1.0e-9)
        @test SciMLBase.successful_retcode(reference)

        @testset "threading = $threading" for threading in (true, PolyesterThreads())
            alg = Alg(threading = threading)
            @test alg.linsolve isa GenericLUFactorization
            sol = solve(prob, alg; abstol = 1.0e-9, reltol = 1.0e-9)
            @test SciMLBase.successful_retcode(sol)
            @test sol.u[end] ≈ reference.u[end] rtol = 1.0e-5
        end
    end
end
