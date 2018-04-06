using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

prob = ODEProblem(DiffEqProblemLibrary.van,[0,2.],(0.0,6),inv(0.003))

# Test if switching back and forth
is_switching_fb(sol) = maximum(diff(find(x->x==2, sol.alg_choice))) > 5
alg = AutoTsit5(Rodas5(); maxstiffstep=5, maxnonstiffstep=5, stiffalgfirst=true)
sol = solve(prob, alg, reltol=1e-5, abstol=1e-5)
@test length(sol.t) < 600
@test typeof(alg.algs[sol.alg_choice[1]]) <: Rodas5
sol = solve(prob, AutoDP5(Rodas5(); maxstiffstep=2, maxnonstiffstep=2,
                          stifftol=(11//10, 9//10), nonstifftol=(9//10, 0//1)),
            reltol=1e-5, abstol=1e-5)
@test length(sol.t) < 610
@test is_switching_fb(sol)
