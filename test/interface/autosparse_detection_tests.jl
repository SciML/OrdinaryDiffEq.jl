using Test, OrdinaryDiffEq, LinearSolve, ADTypes, ForwardDiff, SparseConnectivityTracer, SparseMatrixColorings
import ODEProblemLibrary: prob_ode_2Dlinear
                    

ad = AutoSparse(AutoForwardDiff(), sparsity_detector = TracerSparsityDetector(),
    coloring_algorithm = GreedyColoringAlgorithm())

prob = prob_ode_2Dlinear

@test_nowarn solve(prob, Rodas5(autodiff = ad))

@test_nowarn solve(prob, Rodas5(autodiff = ad, linsolve = LinearSolve.KrylovJL()))

