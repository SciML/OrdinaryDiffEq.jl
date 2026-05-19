using LinearAlgebra, SparseArrays
using OrdinaryDiffEqSDIRK
using SciMLOperators
using RecursiveArrayTools

function rhs_explicit!(du, u, p, t)
    du .= u
end

function rhs_implicit!(du, u, p, t)
    du .= u
end

n = 2
op = MatrixOperator(1.0 * sparse(I, n, n))
# func = SplitFunction(rhs_implicit!, rhs_explicit!)
func = SplitFunction(op, rhs_explicit!)
u0 = ones(n)
u0 = ArrayPartition(u0)
prob = SplitODEProblem{true}(func, u0, (0.0, 1.0))
sol = solve(prob, KenCarp4(), saveat = 0.1)
