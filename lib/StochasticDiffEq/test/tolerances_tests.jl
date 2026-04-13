using StochasticDiffEq, Test, Random
using SDEProblemLibrary: prob_sde_2Dlinear

#=
function f(u,p,t)
  du = similar(u)
  prob_sde_2Dlinear.f(du,u,p,t)
  return du
end
function σ(u,p,t)
  du = similar(u)
  prob_sde_2Dlinear.g(du,u,p,t)
  return du
end

probs = [prob_sde_2Dlinear,
         SDEProblem(f,σ,prob_sde_2Dlinear.u0,prob_sde_2Dlinear.tspan)]
=#

function tolerance_testing(probs, algs)
    for alg in algs, prob in probs

        dt = typeof(alg) <: StochasticDiffEqAdaptiveAlgorithm ? 0.0 : 0.1
        Random.seed!(100)
        sol = solve(prob, alg; dt = dt)

        # Vector of element-wise absolute tolerances
        Random.seed!(100)
        sol2 = solve(prob, alg; dt = dt, abstol = fill(1.0e-2, 4, 2))

        @test sol.t ≈ sol2.t && sol.u ≈ sol2.u

        # Vector of element-wise relative tolerances
        Random.seed!(100)
        sol2 = solve(prob, alg; dt = dt, reltol = fill(1.0e-2, 4, 2))

        @test sol.t ≈ sol2.t && sol.u ≈ sol2.u
    end
    return
end

probs = [prob_sde_2Dlinear]
algs = [SRI(), SRIW1(), SRA1(), SRA(), RKMil(), RKMilGeneral()]
tolerance_testing(probs, algs)
