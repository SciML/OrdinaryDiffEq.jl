using DiffEqBase
using OrdinaryDiffEq

function time_derivative(t,u,du)
  du[1] = -t
end
function time_derivative(::Type{Val{:analytic}},t,u0)
  u0 - t^2/2
end

u0 = [1.0]
tspan = (0.0,1.0)
prob = ODEProblem(time_derivative,u0,tspan)
sol = solve(prob,Rosenbrock32())
sol.errors[:final] < 1e-2
sol = solve(prob,Rosenbrock23())
sol.errors[:final] < 1e-10
sol = solve(prob,ImplicitEuler(),dt=1/10)
sol.errors[:final] < 1e-1
sol = solve(prob,Trapezoid(),dt=1/10)
sol.errors[:final] < 1e-12

for alg in CACHE_TEST_ALGS
  sol = solve(prob,alg,dt=1/100)
  if typeof(alg) <: Euler
    sol.errors[:final] < 1e-3
  else
    sol.errors[:final] < 1e-14
  end
end
