using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

prob = prob_ode_vanderpol

sol = solve(prob, AutoRodas5(alg=Tsit5), dt=1e-5, reltol=1e-5, abstol=1e-5)
@test length(sol.t) < 15
sol = solve(prob, AutoRodas5(alg=DP5), dt=1e-5, reltol=1e-5, abstol=1e-5)
@test length(sol.t) < 15
