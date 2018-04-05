using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

prob = ODEProblem(DiffEqProblemLibrary.van,[0;sqrt(3)],(0.0,3.0),1e6)

# Test if switching back and forth
is_switching_fb(sol) = maximum(diff(find(x->x==2, sol.alg_choice))) > 5
alg = AutoTsit5(Rodas5(); maxstiffstep=15, maxnonstiffstep=8, tol=11//10, stiffalgfirst=true)
sol = solve(prob, alg, dt=1e-5, reltol=1e-5, abstol=1e-5)
@test length(sol.t) < 90
#@test is_switching_fb(sol)
@test typeof(alg.algs[sol.alg_choice[1]]) <: Rodas5
sol = solve(prob, AutoDP5(Rodas5(); maxstiffstep=15, maxnonstiffstep=8, tol=11//10),
            dt=1e-5, reltol=1e-5, abstol=1e-5)
@test length(sol.t) < 150
@test is_switching_fb(sol)
