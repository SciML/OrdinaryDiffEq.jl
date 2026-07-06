using OrdinaryDiffEqSDIRK
using Test

vdp!(du, u, p, t) = (du[1] = u[2]; du[2] = 1.0e3 * (1 - u[1]^2) * u[2] - u[1]; nothing)
prob = ODEProblem(vdp!, [2.0, 0.0], (0.0, 6.3))
ref = solve(prob, ESDIRK659L2SA(predictor = Predictor.Trivial), reltol = 1.0e-12, abstol = 1.0e-12)

methods = [
    ESDIRK325L2SA, ESDIRK54I8L2SA, ESDIRK436L2SA2, ESDIRK437L2SA, ESDIRK547L2SA2,
    ESDIRK659L2SA,
]
predictors = [
    Predictor.Trivial, Predictor.Linear, Predictor.CopyPrev, Predictor.StageExtrap,
    Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder, Predictor.Tableau,
]

@testset "Stage predictors" begin
    for M in methods
        @testset "$(nameof(M))" begin
            triv = solve(prob, M(predictor = Predictor.Trivial), reltol = 1.0e-8, abstol = 1.0e-8)
            for p in predictors
                s = solve(prob, M(predictor = p), reltol = 1.0e-8, abstol = 1.0e-8)
                @test all(isfinite, s.u[end])
                @test maximum(abs.(s.u[end] .- ref.u[end])) < 1.0e-5
            end
            se = solve(prob, M(predictor = Predictor.StageExtrap), reltol = 1.0e-8, abstol = 1.0e-8)
            @test se.stats.nf <= triv.stats.nf
        end
    end
end

# deprecated `extrapolant` keyword still maps onto a predictor
@testset "extrapolant deprecation" begin
    @test ImplicitEuler(extrapolant = :linear).predictor == Predictor.Linear
    @test ImplicitEuler(extrapolant = :constant).predictor == Predictor.Trivial
    @test KenCarp4(extrapolant = :interpolant).predictor == Predictor.MaxOrder
end
