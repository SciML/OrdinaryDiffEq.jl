using OrdinaryDiffEqBDF
using ADTypes
using SparseArrays
import SciMLBase
using Test

# A `DAEFunction` is fully implicit and therefore has no `mass_matrix` field. The sparse
# Jacobian setup shared with the ODE solvers used to read `f.mass_matrix` unconditionally,
# so every fully implicit DAE algorithm errored on a `DAEProblem` carrying a sparse
# `jac_prototype` (SciML/OrdinaryDiffEq.jl#1966).

# Chain of `N` index-1 subsystems:
#   x_i' = -x_i + y_i          (differential)
#   0    =  x_i + y_i - 1      (algebraic)
# with the closed-form solution x_i(t) = (1 - exp(-2t)) / 2, y_i(t) = 1 - x_i(t).
const N = 8

function dae_chain!(res, du, u, p, t)
    for i in 1:N
        res[i] = du[i] + u[i] - u[N + i]
        res[N + i] = u[i] + u[N + i] - 1
    end
    return nothing
end

exact(t) = vcat(fill((1 - exp(-2t)) / 2, N), fill(1 - (1 - exp(-2t)) / 2, N))

const u0 = vcat(zeros(N), ones(N))
const du0 = vcat(ones(N), zeros(N))
const differential_vars = vcat(trues(N), falses(N))
const tspan = (0.0, 1.0)

# Structural pattern of dF/du + γ dF/d(du).
function chain_jac_prototype()
    rows, cols = Int[], Int[]
    for i in 1:N
        append!(rows, (i, i, N + i, N + i))
        append!(cols, (i, N + i, i, N + i))
    end
    return sparse(rows, cols, ones(length(rows)), 2N, 2N)
end

sparse_prob() = DAEProblem(
    DAEFunction(dae_chain!; jac_prototype = chain_jac_prototype()),
    du0, u0, tspan; differential_vars
)
dense_prob() = DAEProblem(dae_chain!, du0, u0, tspan; differential_vars)

# `atol` is per-algorithm because the family spans orders: `DImplicitEuler` is first
# order, so it cannot reach the accuracy the adaptive-order `DFBDF` does on this span.
const CASES = (
    ("DFBDF", DFBDF(), 1.0e-6),
    ("DFBDF, AutoFiniteDiff", DFBDF(autodiff = AutoFiniteDiff()), 1.0e-6),
    ("DImplicitEuler", DImplicitEuler(), 1.0e-5),
    ("DABDF2", DABDF2(), 1.0e-6),
)

@testset "$name with sparse jac_prototype" for (name, alg, atol) in CASES
    sparse_sol = solve(sparse_prob(), alg; abstol = 1.0e-10, reltol = 1.0e-10)
    @test SciMLBase.successful_retcode(sparse_sol)
    @test maximum(maximum(abs, s .- exact(t)) for (s, t) in zip(sparse_sol.u, sparse_sol.t)) < atol

    # Declaring the sparsity must not change the trajectory the dense path produces.
    dense_sol = solve(dense_prob(), alg; abstol = 1.0e-10, reltol = 1.0e-10)
    @test sparse_sol(tspan[2]) ≈ dense_sol(tspan[2]) rtol = 1.0e-8
end
