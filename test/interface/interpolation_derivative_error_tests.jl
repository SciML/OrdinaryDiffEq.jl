using OrdinaryDiffEq
function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan)

out = rand(3)

for alg in [Rosenbrock23(), Rodas4(), Rodas5P(), Tsit5(), DP5(),
    BS5(), OwrenZen3(), OwrenZen4(), OwrenZen5(),
    Vern6(), Vern7(), Vern8(), Vern9(), BS3()]
    sol = solve(prob, alg)
    @test_throws OrdinaryDiffEq.DerivativeOrderNotPossibleError sol(0.5, Val{10})
    @test_throws OrdinaryDiffEq.DerivativeOrderNotPossibleError sol(out, 0.5, Val{10})
end
