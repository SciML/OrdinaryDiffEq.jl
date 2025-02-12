using OrdinaryDiffEq, LinearSolve, ADTypes, ForwardDiff, SparseConnectivityTracer, SparseMatrixColorings
import ODEProblemLibrary: prob_ode_2Dlinear
                    

ad = AutoSparse(AutoForwardDiff(), sparsity_detector = TracerSparsityDetector(),
    coloring_algorithm = GreedyColoringAlgorithm())

prob = prob_ode_2Dlinear

@test_no_warn solve(prob, Rodas5(autodiff = ad))

@test_no_warn solve(prob, Rodas5(autodiff = ad, linsolve = LinearSolve.KrylovJL()))

