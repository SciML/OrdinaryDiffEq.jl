using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, Test
using SciMLBase: DAEProblem

function dae_lin!(res, du, u, p, t)
    @. res = du - u
    return nothing
end
dae_lin(du, u, p, t) = @. du - u

const prob_iip = DAEProblem(dae_lin!, [1.0, 1.0], [1.0, 1.0], (0.0, 1.0))
const prob_oop = DAEProblem{false}(dae_lin, [1.0, 1.0], [1.0, 1.0], (0.0, 1.0))
const dae_bdf_algs = (DImplicitEuler(), DABDF2(), DFBDF(), DNordsieckBDF())

@testset "DAE BDF uses COEFFICIENT_MULTISTEP (issue #4113)" begin
    for (prob, label) in ((prob_iip, "iip"), (prob_oop, "oop"))
        @testset "$label" begin
            for alg in dae_bdf_algs
                integrator = init(
                    prob, alg; dt = 0.1, adaptive = false, save_everystep = false
                )
                step!(integrator)
                step!(integrator)
                nlsolver = integrator.cache.nlsolver
                @test nlsolver.method === OrdinaryDiffEqCore.COEFFICIENT_MULTISTEP
                @test nlsolver.z ≈ integrator.u
            end
        end
    end
end
