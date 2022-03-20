# destats.nf tests
using OrdinaryDiffEq, Test
x = Ref(0)
function f(u,p,t)
  x[] += 1
  return 5*u
end
u0 = [1.0, 1.0]
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)

x[] = 0
sol = solve(prob, Vern7())
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, Vern8())
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, Vern9())
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, Tsit5())
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, BS3())
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, KenCarp4(;autodiff=true))
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, KenCarp4(;autodiff=false, diff_type=Val{:forward}))
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, KenCarp4(;autodiff=false, diff_type=Val{:central}))
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, KenCarp4(;autodiff=false, diff_type=Val{:complex}))
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, Rosenbrock23(;autodiff=true))
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, Rosenbrock23(;autodiff=false, diff_type=Val{:forward}))
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, Rosenbrock23(;autodiff=false, diff_type=Val{:central}))
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, Rosenbrock23(;autodiff=false, diff_type=Val{:complex}))
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, Rodas5(;autodiff=true))
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, Rodas5(;autodiff=false, diff_type=Val{:forward}))
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, Rodas5(;autodiff=false, diff_type=Val{:central}))
@test x[] == sol.destats.nf

x[] = 0
sol = solve(prob, Rodas5(;autodiff=false, diff_type=Val{:complex}))
@test x[] == sol.destats.nf

function g(du, u,p,t)
  x[] += 1
  @. du = 5*u
end
probip = ODEProblem(g,u0,tspan)

x[] = 0
sol = solve(probip, ROCK4())
@test x[] == sol.destats.nf
