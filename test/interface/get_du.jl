using OrdinaryDiffEq, OrdinaryDiffEqCore, Test
function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 3.0)
condition(u, t, integrator) = t == 0.2
cache = zeros(3)
affect!(integrator) = cache .= get_du(integrator)
dusave = DiscreteCallback(condition, affect!)
affect2!(integrator) = get_du!(cache, integrator)
dusave_inplace = DiscreteCallback(condition, affect2!)

prob = ODEProblem(lorenz!, u0, tspan)
sol = solve(
    prob, Tsit5(), tstops = [0.2], callback = dusave, abstol = 1.0e-12, reltol = 1.0e-12
)
res = copy(cache)

for alg in [
        Vern6(), Vern7(), Vern8(), Vern9(), Rodas4(), Rodas4P(),
        Rodas5(), Rodas5P(), TRBDF2(), KenCarp4(), FBDF(), QNDF(),
    ]
    sol = solve(
        prob, alg, tstops = [0.2], callback = dusave, abstol = 1.0e-12, reltol = 1.0e-12
    )
    @test res ≈ cache rtol = 1.0e-5

    sol = solve(
        prob, alg, tstops = [0.2], callback = dusave_inplace,
        abstol = 1.0e-12, reltol = 1.0e-12
    )
    @test res ≈ cache rtol = 1.0e-5
end
