using OrdinaryDiffEq, Test

f_2dlinear = (du, u, p, t) -> (@. du = p * u)

prob_ode_2Dlinear = ODEProblem(f_2dlinear, rand(4, 2), (0.0, 1.0), 1.01)
sol = solve(prob_ode_2Dlinear)

tsitsol = solve(prob_ode_2Dlinear, Tsit5())
# test that default isn't much worse than Tsit5 (we expect it to use Tsit5 for this).
@test sol.stats.naccept < tsitsol.stats.naccept + 2
@test sol.stats.nf < tsitsol.stats.nf + 20

sol = solve(prob_ode_2Dlinear, reltol=1e-10)
vernsol = solve(prob_ode_2Dlinear, Vern7(), reltol=1e-10)
# test that default isn't much worse than Tsit5 (we expect it to use Tsit5 for this).
@test sol.stats.naccept < tsitsol.stats.naccept + 2
@test sol.stats.nf < tsitsol.stats.nf + 20

function rober(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    [-k₁ * y₁ + k₃ * y₂ * y₃,
      k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2,
      k₂ * y₂^2]
end

prob_rober = ODEProblem(rober, [1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
sol = solve(prob_rober)
rosensol = solve(prob_rober, Rosenbrock23())
