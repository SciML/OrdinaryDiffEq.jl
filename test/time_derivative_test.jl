using OrdinaryDiffEq, Base.Test

function time_derivative(t,u,du)
  du[1] = -t
end
function time_derivative(::Type{Val{:analytic}},t,u0)
  u0 - t^2/2
end

u0 = [1.0]
tspan = (0.0,1.0)
prob = ODEProblem(time_derivative,u0,tspan)
sol = solve(prob,Rosenbrock32(),reltol=1e-9,abstol=1e-9)
@test sol.errors[:final] < 1e-5
sol = solve(prob,Rosenbrock23())
@test sol.errors[:final] < 1e-10
sol = solve(prob,GenericImplicitEuler(),dt=1/10)
@test sol.errors[:final] < 1e-1
sol = solve(prob,GenericTrapezoid(),dt=1/10)
@test sol.errors[:final] < 1e-12
sol = solve(prob,Euler(),dt=1/100)
@test sol.errors[:final] < 6e-3

for alg in CACHE_TEST_ALGS
  sol = solve(prob,alg,dt=1/10)
  if !(typeof(alg) <: Euler)
    @test sol.errors[:final] < 4e-14
  end
end
