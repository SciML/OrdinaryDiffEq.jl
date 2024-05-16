using OrdinaryDiffEq, Test, LinearSolve

f_2dlinear = (du, u, p, t) -> (@. du = p * u)

prob_ode_2Dlinear = ODEProblem(f_2dlinear, rand(4, 2), (0.0, 1.0), 1.01)
sol = solve(prob_ode_2Dlinear)

tsitsol = solve(prob_ode_2Dlinear, Tsit5())
# test that default isn't much worse than Tsit5 (we expect it to use Tsit5 for this).
@test sol.stats.naccept == tsitsol.stats.naccept
@test sol.stats.nf == tsitsol.stats.nf
@test all(isequal(1), sol.alg_choice)

sol = solve(prob_ode_2Dlinear, reltol=1e-10)
vernsol = solve(prob_ode_2Dlinear, Vern7(), reltol=1e-10)
# test that default isn't much worse than Tsit5 (we expect it to use Tsit5 for this).
@test sol.stats.naccept == vernsol.stats.naccept
@test sol.stats.nf == vernsol.stats.nf
@test all(isequal(2), sol.alg_choice)

function rober(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    [-k₁ * y₁ + k₃ * y₂ * y₃,
      k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2,
      k₂ * y₂^2]
end
prob_rober = ODEProblem(rober, [1.0,0.0,0.0],(0.0,1e3),(0.04,3e7,1e4))
sol = solve(prob_rober)
rosensol = solve(prob_rober, AutoTsit5(Rosenbrock23()))
# test that default has the same performance as AutoTsit5(Rosenbrock23()) (which we expect it to use for this).
@test sol.stats.naccept == rosensol.stats.naccept
@test sol.stats.nf == rosensol.stats.nf
@test unique(sol.alg_choice) == [1,3]
@test sol.alg_choice[1] == 1
@test sol.alg_choice[end] == 3

sol = solve(prob_rober, reltol=1e-7, abstol=1e-7)
rosensol = solve(prob_rober, AutoVern7(Rodas5P()), reltol=1e-7, abstol=1e-7)
# test that default has the same performance as AutoTsit5(Rosenbrock23()) (which we expect it to use for this).
@test sol.stats.naccept == rosensol.stats.naccept
@test sol.stats.nf == rosensol.stats.nf
@test unique(sol.alg_choice) == [2,4]
@test sol.alg_choice[1] == 2
@test sol.alg_choice[end] == 4

function exrober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du .= vcat([-k₁ * y₁ + k₃ * y₂ * y₃,
      k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2,
      k₂ * y₂^2, ], u[4:end])
end

for n in (100, 1000)
    stiffalg = n < 500 ? 5 : 6
    linsolve = n < 500 ? nothing : KrylovJL_GMRES()
    prob_ex_rober = ODEProblem(exrober, vcat([1.0,0.0,0.0], ones(n)),(0.0,100.0),(0.04,3e7,1e4))
    sol = solve(prob_ex_rober)
    fsol = solve(prob_ex_rober, AutoTsit5(FBDF(;linsolve)))
    # test that default has the same performance as AutoTsit5(Rosenbrock23()) (which we expect it to use for this).
    @test sol.stats.naccept == fsol.stats.naccept
    @test sol.stats.nf == fsol.stats.nf
    @test unique(sol.alg_choice) == [1,stiffalg]
    @test sol.alg_choice[1] == 1
    @test sol.alg_choice[end] == stiffalg
end
