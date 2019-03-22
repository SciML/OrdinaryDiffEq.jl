using OrdinaryDiffEq, Test
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: van

prob1 = ODEProblem(van,  [0,2.],(0.0,6),inv(0.003))
prob2 = ODEProblem(van.f,[0,2.],(0.0,6),inv(0.003))
# out-of-place test
function _van(u, p, t)
  μ = p[1]
  [μ*((1-u[2]^2)*u[1] - u[2]),
   1*u[1]]
end
prob3 = ODEProblem(_van,[0,2.],(0.0,6),inv(0.003))
probArr = [prob1, prob2, prob3]
# Test if switching back and forth
is_switching_fb(sol) = maximum(diff(findall(x->x==2, sol.alg_choice))) > 5
for prob in probArr
  sol = @test_nowarn solve(prob, AutoTsit5(Rosenbrock23(autodiff=false)))
  @test is_switching_fb(sol)
  alg = AutoTsit5(Rodas5(); maxstiffstep=5, maxnonstiffstep=5, stiffalgfirst=true)
  sol = solve(prob, alg)
  @test length(sol.t) < 280
  @test typeof(alg.algs[sol.alg_choice[1]]) <: Rodas5
  @test is_switching_fb(sol)
  sol = solve(prob, AutoDP5(Rodas5(); maxstiffstep=2, maxnonstiffstep=2,
                            stifftol=11//10, nonstifftol=9//10),
                            reltol=1e-5, abstol=1e-5)
  @test length(sol.t) < 625
  @test is_switching_fb(sol)

  sol = solve(prob,AutoVern6(Kvaerno3(); maxstiffstep=4, maxnonstiffstep=4))
  @test length(sol.t) < 690
  @test is_switching_fb(sol)
  sol = solve(prob,AutoVern7(Hairer42(); maxstiffstep=4, maxnonstiffstep=4))
  @test length(sol.t) < 540
  @test is_switching_fb(sol)
  sol = solve(prob,AutoVern8(Rosenbrock23(); maxstiffstep=4, maxnonstiffstep=4))
  @test length(sol.t) < 910
  @test is_switching_fb(sol)
  sol = solve(prob,AutoVern9(KenCarp3(); maxstiffstep=4, maxnonstiffstep=4))
  @test length(sol.t) < 470
  @test is_switching_fb(sol)
  sol = solve(prob,AutoVern9(KenCarp3(autodiff=false); maxstiffstep=4, maxnonstiffstep=4))
  @test length(sol.t) < 470
  @test is_switching_fb(sol)
end
