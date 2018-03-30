using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

prob = DiffEqProblemLibrary.prob_ode_vanstiff

sol = solve(prob, AutoRodas5(Tsit5(); maxstiffstep=15), dt=1e-5, reltol=1e-5, abstol=1e-5)
@test length(sol.t) < 40
@test length(unique(sol.alg_choice)) == 2
sol = solve(prob, AutoRodas5(DP5(); maxstiffstep=15), dt=1e-5, reltol=1e-5, abstol=1e-5)
@test length(sol.t) < 42
@test length(unique(sol.alg_choice)) == 2
