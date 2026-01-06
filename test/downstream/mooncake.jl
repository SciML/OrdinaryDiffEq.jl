using Mooncake, OrdinaryDiffEq, StaticArrays

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

const _saveat = SA[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]

function f(u0::Array{Float64})
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = SciMLBase.solve(prob, Tsit5(), saveat = _saveat, sensealg = SciMLBase.SensitivityADPassThrough())
    return sum(sol)
end;
u0 = [1.0; 0.0; 0.0]
mooncake_gradient(f, x) = Mooncake.value_and_gradient!!(Mooncake.build_rrule(f, x), f, x)[2][2]
@test_broken mooncake_gradient(f, u0)
