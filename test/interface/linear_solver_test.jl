using Test, OrdinaryDiffEq

using SparseArrays
using LinearAlgebra, Random
N = 30
AA = sprand(MersenneTwister(12),  N, N, 0.5)
mm = sprand(MersenneTwister(123), N, N, 0.5)
A = DiffEqArrayOperator(AA)
M = DiffEqArrayOperator(mm'mm)
u0 = ones(N)
prob = ODEProblem(ODEFunction(A; mass_matrix=M), u0, (0.0, 1.0))

for alg in [Rosenbrock23()]
  sol = solve(prob, alg)
  @test sol.destats.njacs == 0
  @test sol.destats.nw == 1
end
