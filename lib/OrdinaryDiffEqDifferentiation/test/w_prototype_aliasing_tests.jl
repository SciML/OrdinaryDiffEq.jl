using OrdinaryDiffEqSDIRK, OrdinaryDiffEqDifferentiation, SciMLOperators, LinearAlgebra
using SciMLBase, Test

@testset "integrators from one ODEFunction do not share the W_prototype" begin
    mutable struct ParamA
        a::Float64
    end
    f!(du, u, p, t) = (du .= p.a .* u)
    Jop = MatrixOperator(
        zeros(2, 2); update_func! = (A, u, p, t) -> (A .= p.a * Matrix(1.0I, 2, 2))
    )
    Wop = MatrixOperator(
        zeros(2, 2);
        update_func! = (A, u, p, t; gamma = 1.0) ->
        (A .= p.a * Matrix(1.0I, 2, 2) .- Matrix(1.0I, 2, 2) ./ gamma),
        accepted_kwargs = (:gamma,),
    )
    odef = ODEFunction(f!; jac_prototype = Jop, W_prototype = Wop)
    prob1 = ODEProblem(odef, ones(2), (0.0, 1.0), ParamA(-1.0))
    prob2 = ODEProblem(odef, ones(2), (0.0, 1.0), ParamA(-1000.0))

    i1 = init(prob1, TRBDF2(); abstol = 1.0e-8, reltol = 1.0e-8)
    i2 = init(prob2, TRBDF2(); abstol = 1.0e-8, reltol = 1.0e-8)
    @test i1.cache.nlsolver.cache.W !== i2.cache.nlsolver.cache.W
    @test i1.cache.nlsolver.cache.W !== odef.W_prototype
    @test i1.cache.nlsolver.cache.J !== odef.jac_prototype

    for _ in 1:200
        i1.t < 1.0 && step!(i1)
        i2.t < 1.0 && step!(i2)
        i1.t >= 1.0 && i2.t >= 1.0 && break
    end
    @test i1.u[1] ≈ exp(-1.0) rtol = 1.0e-3
    @test abs(i2.u[1]) < 1.0e-6
end
