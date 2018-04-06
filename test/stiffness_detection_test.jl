using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

prob = ODEProblem(DiffEqProblemLibrary.van,[0,2.],(0.0,6),inv(0.003))

# Test if switching back and forth
is_switching_fb(sol) = maximum(diff(find(x->x==2, sol.alg_choice))) > 5
alg = AutoTsit5(Rodas5(); maxstiffstep=5, maxnonstiffstep=5, stiffalgfirst=true)
sol = solve(prob, alg)
@test length(sol.t) < 400
@test typeof(alg.algs[sol.alg_choice[1]]) <: Rodas5
@test is_switching_fb(sol)
sol = solve(prob, AutoDP5(Rodas5(); maxstiffstep=2, maxnonstiffstep=2,
                          stifftol=11//10, nonstifftol=9//10),
                          reltol=1e-5, abstol=1e-5)
@test length(sol.t) < 610
@test is_switching_fb(sol)

sol = solve(prob,AutoVern6(Rosenbrock23(); maxstiffstep=4, maxnonstiffstep=4))
@test length(sol.t) < 600
@test is_switching_fb(sol)
sol = solve(prob,AutoVern7(Rosenbrock23(); maxstiffstep=4, maxnonstiffstep=4))
@test length(sol.t) < 600
@test is_switching_fb(sol)
sol = solve(prob,AutoVern8(Rosenbrock23(); maxstiffstep=4, maxnonstiffstep=4))
@test length(sol.t) < 600
@test is_switching_fb(sol)
sol = solve(prob,AutoVern9(Rosenbrock23(); maxstiffstep=4, maxnonstiffstep=4))
@test length(sol.t) < 600
@test is_switching_fb(sol)
