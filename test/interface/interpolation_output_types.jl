using OrdinaryDiffEq, RecursiveArrayTools, Test

# in terms of the voltage across all three elements
rlc1!(v′, v, (R, L, C), t) = -(v′ / R + v / L) / C
identity_f(v, u, p, t) = v # needed to form second order dynamical ODE

function setup_rlc(R, L, C; v_init = 0.0, v′_init = 0.0, tspan = (0.0, 50.0))
    DynamicalODEProblem{false}(rlc1!, identity_f, v′_init, v_init, tspan, (R, L, C))
end

# simulate voltage impulse
R, L, C = 10, 0.3, 2

prob = setup_rlc(R, L, C, v_init = 2.0)

res1 = solve(prob, Vern8(), dt = 1 / 10, saveat = 1 / 10)
res3 = solve(prob, CalvoSanz4(), dt = 1 / 10, saveat = 1 / 10)

sol = solve(prob, CalvoSanz4(), dt = 1 / 10)
@test sol(0.32) isa RecursiveArrayTools.ArrayPartition
@test sol(0.32, Val{1}) isa RecursiveArrayTools.ArrayPartition
@test sol(0.32, Val{2}) isa RecursiveArrayTools.ArrayPartition
@test sol(0.32, Val{3}) isa RecursiveArrayTools.ArrayPartition

function f(du, u, p, t)
    du .= u
    nothing
end
dprob = DiscreteProblem(f, [1, 2, 3], (0, 100))
sol = solve(dprob, FunctionMap())
@test sol(0:0.1:100; idxs = 1) isa DiffEqArray
@test length(sol(0:0.1:100; idxs = 1)) == length(0:0.1:100)
@test length(sol(0:0.1:100; idxs = 1).u[1]) == 1
sol(0:0.1:100; idxs = [1, 2])

@test sol(0:0.1:100; idxs = [1, 2]) isa DiffEqArray
@test length(sol(0:0.1:100; idxs = [1, 2])) == length(0:0.1:100)
@test length(sol(0:0.1:100; idxs = [1, 2]).u[1]) == 2

@test sol(0:0.1:100) isa DiffEqArray
@test length(sol(0:0.1:100)) == length(0:0.1:100)
@test length(sol(0:0.1:100).u[1]) == 3

## Test DPRKN Interpolation

#Parameters
ω = 1

#Initial Conditions
x₀ = [0.0]
dx₀ = [π / 2]
tspan = (0.0, 2π)

ϕ = atan((dx₀[1] / ω) / x₀[1])
A = √(x₀[1]^2 + dx₀[1]^2)

function harmonicoscillator(ddu, du, u, ω, t)
    ddu .= -ω^2 * u
end

prob = SecondOrderODEProblem(harmonicoscillator, dx₀, x₀, tspan, ω)
sol = solve(prob, DPRKN6())
@test sol(0.5) isa RecursiveArrayTools.ArrayPartition
@test sol(0.5; idxs = 1) isa Number
@test sol(0.5; idxs = [1]) isa Vector
@test sol(0.5; idxs = [1, 2]) isa Vector
@test Vector(sol(0.5)) == sol(0.5; idxs = [1, 2])
@test reverse(Vector(sol(0.5))) == sol(0.5; idxs = [2, 1])
@test Vector(sol(0.5)) == [sol(0.5; idxs = 1); sol(0.5; idxs = 2)]
