using OrdinaryDiffEq, Test
function some_arbitrary_function!(du, u, p, τ)
    du = u / 100
    return nothing #function returns nothing
end
NF = Float32                   # Number format
u = NF[0.0, 100.0, π / 2.0, 0.0] #Initial conditions e.g. t, r θ, ϕ
params = NF[0.1, 3.5f-7]  #Some arbitrary parameters, unused in this example, just a placeholder
tspan = (zero(NF), NF(1.0e3)) #integrate from t=0 to t = 1000

ode_prob = ODEProblem(some_arbitrary_function!, u, tspan, params)

for alg in [
        Euler(), Midpoint(), Heun(), Ralston(), RK4(), SSPRK104(), SSPRK22(), SSPRK33(),
        SSPRK43(), SSPRK432(), BS3(), BS5(), DP5(), DP8(), Feagin10(), Feagin12(),
        Feagin14(), TanYam7(), Tsit5(), TsitPap8(), Vern6(), Vern7(), Vern8(), Vern9(),
    ]
    @test ode_solution = solve(ode_prob, alg, dt = 1.0f-1).retcode === ReturnCode.Success
end
