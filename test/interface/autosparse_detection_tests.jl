using Test, OrdinaryDiffEq, LinearSolve, ADTypes, ForwardDiff, SparseConnectivityTracer,
    SparseMatrixColorings
import ODEProblemLibrary: prob_ode_2Dlinear

ad = AutoSparse(
    AutoForwardDiff(), sparsity_detector = TracerSparsityDetector(),
    coloring_algorithm = GreedyColoringAlgorithm()
)

prob = prob_ode_2Dlinear

@test_nowarn solve(prob, Rodas5P(autodiff = ad))

@test_nowarn solve(prob, Rodas5P(autodiff = ad, linsolve = LinearSolve.KrylovJL_GMRES()))

@test_nowarn solve(prob, FBDF(autodiff = ad))

# Test that no dense matrices are made sparse
diag_prob = ODEProblem((du, u, p, t) -> du .= -1.0 .* u, rand(Int(1.0e7)), (0, 1.0))

@test_nowarn solve(diag_prob, Rodas5P(autodiff = ad, linsolve = LinearSolve.KrylovJL_GMRES()))
