using StochasticDiffEq, Test, Random
using SDEProblemLibrary: prob_sde_2Dlinear
Random.seed!(100)
prob = prob_sde_2Dlinear

## Solve and plot
sol1 = solve(prob, RKMil(), abstol = 1, reltol = 0, adaptive = true)
err1 = sol1.errors[:final]

sol1 = solve(prob, RKMil(), abstol = 1.0e-3, reltol = 1.0e-3, adaptive = true)
err12 = sol1.errors[:final]
@test err12 < err1

sol1 = solve(prob, RKMilGeneral(), abstol = 1, reltol = 0, adaptive = true)
err11 = sol1.errors[:final]

sol1 = solve(prob, RKMilGeneral(), adaptive = true)
err112 = sol1.errors[:final]
@test err112 < err11

sol2 = solve(prob, SRI(), dt = 1 / 2^(4), abstol = 1, reltol = 0)
err2 = sol2.errors[:final]

println("1e-1")
sol3 = solve(prob, SRI(), dt = 1 / 2^(4), abstol = 1.0e-1, reltol = 0)
err3 = sol3.errors[:final]

println("1e-2")
sol4 = solve(prob, SRI(), dt = 1 / 2^(4), abstol = 1.0e-2, reltol = 0)
err4 = sol4.errors[:final]
@test err4 < err2
println("1e-3")
sol5 = solve(prob, SRI(), dt = 1 / 2^(4), abstol = 1.0e-3, reltol = 0)
err5 = sol5.errors[:final]
@test err5 < err3
println(
    """
    Final error for the solutions were:
              $err1
              $err11
              $err2
              $err3
              $err4
              $err5
              """
)

sol6 = solve(prob, SRI())
