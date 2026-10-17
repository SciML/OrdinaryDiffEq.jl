using OrdinaryDiffEqBDF, LinearSolve, SparseArrays, Test
using LinearAlgebra: I, lu, mul!, norm
using ADTypes: AutoForwardDiff, AutoFiniteDiff, AutoSparse
using SciMLBase: DAEProblem, successful_retcode
using SparseConnectivityTracer: TracerSparsityDetector
using SparseMatrixColorings: GreedyColoringAlgorithm

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/1673

function decay!(res, du, u, p, t)
    @. res = du + u
    return nothing
end

# Semi-explicit index-1 DAE: u[3] is algebraic, so ∂F/∂du is singular.
function semi_explicit!(res, du, u, p, t)
    res[1] = du[1] + u[1] - u[2]^2
    res[2] = du[2] + 2u[2] - u[3]
    res[3] = u[3] - sin(u[1])
    return nothing
end

@testset "DFBDF with a matrix-free linear solver on a DAEProblem" begin
    prob = DAEProblem(decay!, -ones(4), ones(4), (0.0, 1.0))
    sol = solve(prob, DFBDF(; linsolve = KrylovJL_GMRES()); abstol = 1.0e-8, reltol = 1.0e-8)
    @test successful_retcode(sol)
    @test sol.u[end] ≈ fill(exp(-1), 4) rtol = 1.0e-5
    # Matrix-free: no Jacobian is ever formed.
    @test sol.stats.njacs == 0

    u0 = [1.0, 0.5, sin(1.0)]
    du0 = [-1.0 + 0.25, -1.0 + sin(1.0), 0.0]
    prob = DAEProblem(
        semi_explicit!, du0, u0, (0.0, 2.0); differential_vars = [true, true, false]
    )
    opts = (; abstol = 1.0e-9, reltol = 1.0e-9)
    reference = solve(prob, DFBDF(); opts...)
    for autodiff in (AutoForwardDiff(), AutoFiniteDiff())
        sol = solve(prob, DFBDF(; autodiff, linsolve = KrylovJL_GMRES()); opts...)
        @test successful_retcode(sol)
        @test sol.stats.njacs == 0
        @test sol.u[end] ≈ reference.u[end] rtol = 1.0e-6
    end

    # User `precs` are wrapped in the residual weights, so the weighted Krylov `abstol`
    # holds for them as it does for the default preconditioner.
    linsolve = KrylovJL_GMRES(; precs = Returns((I, I)))
    sol = solve(prob, DFBDF(; linsolve); opts...)
    @test successful_retcode(sol)
    @test sol.u[end] ≈ reference.u[end] rtol = 1.0e-6
end

# A linear index-1 DAE whose residual rows are rotated against the unknowns, the way a
# ModelingToolkit discretization can order them: every diagonal entry of
# ∂F/∂u + cj * ∂F/∂(du) is zero.
function rotated_chain!(res, du, u, p, t)
    n = length(u)
    for i in 1:(n - 1)
        res[i + 1] = du[i] + 2u[i] - u[i + 1]
    end
    res[1] = u[n] - u[1] / 2
    return nothing
end

@testset "DFBDF with a preconditioned Krylov solve on a DAEProblem" begin
    n = 40
    u0 = [ones(n - 1); 0.5]
    du0 = [fill(-1.0, n - 2); -1.5; -0.5]
    prob = DAEProblem(
        rotated_chain!, du0, u0, (0.0, 1.0); differential_vars = [trues(n - 1); false]
    )
    opts = (; abstol = 1.0e-9, reltol = 1.0e-9)
    reference = solve(prob, DFBDF(); opts...)
    sparse_fd = AutoSparse(
        AutoForwardDiff(); sparsity_detector = TracerSparsityDetector(),
        coloring_algorithm = GreedyColoringAlgorithm()
    )
    for (autodiff, matrix_type) in ((AutoForwardDiff(), Matrix), (sparse_fd, SparseMatrixCSC))
        calls = Ref(0)
        mismatch = Ref(0.0)
        iterations = Ref(0)
        # `convert(AbstractMatrix, W)` is the assembled `J_u + cj * J_du` that the
        # operator applies matrix-free. The first call, from `init`, precedes every
        # Jacobian evaluation and sees placeholder values.
        precs = function (W, p)
            calls[] += 1
            Wc = convert(AbstractMatrix, W)
            @test Wc isa matrix_type
            if calls[] > 1
                v = collect(1.0:n)
                mismatch[] = max(mismatch[], norm(mul!(similar(v), W, v) - Wc * v) / norm(v))
            end
            return (lu(Wc), I)
        end
        linsolve = KrylovJL_GMRES(; precs, callback = _ -> (iterations[] += 1; false))
        sol = solve(prob, DFBDF(; autodiff, concrete_jac = true, linsolve); opts...)
        @test successful_retcode(sol)
        @test sol.u[end] ≈ reference.u[end] rtol = 1.0e-7
        @test sol.stats.njacs > 0
        @test calls[] > 1
        @test mismatch[] < 1.0e-12
        # An exact preconditioner leaves GMRES about one iteration per solve.
        @test iterations[] / sol.stats.nsolve < 2
    end
end
