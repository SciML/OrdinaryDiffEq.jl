using OrdinaryDiffEq, Test

function time_derivative(du,u,p,t)
  du[1] = -t
end
function time_derivative_analytic(u0,p,t)
  u0 .- t.^2 ./ 2
end

ff_time_derivative = ODEFunction(time_derivative,
                                 analytic=time_derivative_analytic)
u0 = [1.0]
tspan = (0.0,1.0)
prob = ODEProblem(ff_time_derivative,u0,tspan)

sol = solve(prob,Rosenbrock32(),reltol=1e-9,abstol=1e-9)
@test sol.errors[:final] < 1e-5
@show Rosenbrock23()
sol = solve(prob,Rosenbrock23(),reltol=1e-9,abstol=1e-9)
@test sol.errors[:final] < 1e-10
@show Rodas4
sol = solve(prob,Rodas4(),reltol=1e-9,abstol=1e-9)
@test sol.errors[:final] < 1e-10
@show Rodas5
sol = solve(prob,Rodas5(),reltol=1e-9,abstol=1e-9)
@test sol.errors[:final] < 1e-10
@show Veldd4
sol = solve(prob,Veldd4(),reltol=1e-9,abstol=1e-9)
@test sol.errors[:final] < 1e-10
@show KenCarp3
sol = solve(prob,KenCarp3(),reltol=1e-12,abstol=1e-12)
@test length(sol) > 2
@test sol.errors[:final] < 1e-10
@show KenCarp4
sol = solve(prob,KenCarp4(),reltol=1e-12,abstol=1e-12)
@test length(sol) > 2
@test sol.errors[:final] < 1e-10
@show KenCarp5
sol = solve(prob,KenCarp5(),reltol=1e-12,abstol=1e-12)
@test length(sol) > 2
@test sol.errors[:final] < 1e-10
@show TRBDF2
sol = solve(prob,TRBDF2(),reltol=1e-9,abstol=1e-9)
@test sol.errors[:final] < 1e-10
sol = solve(prob,ImplicitEuler(),dt=1/10)
@test sol.errors[:final] < 1e-1
sol = solve(prob,Trapezoid(),dt=1/10)
@test sol.errors[:final] < 1e-12
sol = solve(prob,Euler(),dt=1/100)
@test sol.errors[:final] < 6e-3


const CACHE_TEST_ALGS = [Euler(),Midpoint(),RK4(),SSPRK22(),SSPRK33(),SSPRK53(),
  SSPRK63(),SSPRK73(),SSPRK83(),SSPRK432(),SSPRK932(),SSPRK54(),SSPRK104(),CarpenterKennedy2N54(),
  BS3(),BS5(),DP5(),DP8(),Feagin10(),Feagin12(),Feagin14(),TanYam7(),
  Tsit5(),TsitPap8(),Vern6(),Vern7(),Vern8(),Vern9(),OwrenZen3(),OwrenZen4(),OwrenZen5()]

for alg in CACHE_TEST_ALGS
  sol = solve(prob,alg,dt=1/10)
  if !(typeof(alg) <: Euler)
    @test sol.errors[:final] < 4e-14
  end
end
