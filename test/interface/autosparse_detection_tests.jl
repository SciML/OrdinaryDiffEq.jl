using Test, OrdinaryDiffEq, LinearSolve, ADTypes, ForwardDiff, SparseConnectivityTracer,
      SparseMatrixColorings
import ODEProblemLibrary: prob_ode_2Dlinear

ad = AutoSparse(AutoForwardDiff(), sparsity_detector = TracerSparsityDetector(),
    coloring_algorithm = GreedyColoringAlgorithm())

prob = prob_ode_2Dlinear

@test_nowarn solve(prob, Rodas5P(autodiff = ad))

@test_nowarn solve(prob, Rodas5P(autodiff = ad, linsolve = LinearSolve.KrylovJL_GMRES()))

# Test that no dense matrices are made sparse
diag_prob = ODEProblem((du, u, p, t) -> du .= -1.0 .* u, rand(Int(1e7)), (0, 1.0))

solve(diag_prob, Rodas5P(autodiff = ad), tspan = (0.0, 1e-8))

@test_nowarn solve(diag_prob, Rodas5P(autodiff = ad))

@test_nowarn solve(
    diag_prob, Rodas5P(auodiff = ad, linsolve = LinearSolve.KrylovJL_GMRES()))