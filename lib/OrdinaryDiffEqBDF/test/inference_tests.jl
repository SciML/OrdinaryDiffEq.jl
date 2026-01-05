using OrdinaryDiffEqBDF, ADTypes, Test
using NonlinearSolve: TrustRegion

prob = ODEProblem((du, u, p, t) -> du .= u, zeros(1), (0.0, 1.0))
nlalg = FBDF(
    autodiff = false,
    nlsolve = OrdinaryDiffEqBDF.NonlinearSolveAlg(TrustRegion(autodiff = AutoFiniteDiff()))
)
basicalg = FBDF(autodiff = false)
basicalgad = FBDF()

nlsolver = @inferred OrdinaryDiffEqBDF.build_nlsolver(
    basicalg, prob.u0, prob.u0, prob.p, 0.0, 0.0, prob.f, prob.u0, Float64,
    Float64, Float64, 0.0, 0.0, Val(true)
)
nlsolver = @inferred OrdinaryDiffEqBDF.build_nlsolver(
    nlalg, prob.u0, prob.u0, prob.p, 0.0, 0.0, prob.f, prob.u0, Float64,
    Float64, Float64, 0.0, 0.0, Val(true)
)
nlsolver = @test_throws Any @inferred OrdinaryDiffEqBDF.build_nlsolver(
    basicalgad, prob.u0, prob.u0, prob.p, 0.0, 0.0, prob.f, prob.u0, Float64,
    Float64, Float64, 0.0, 0.0, Val(true)
)
