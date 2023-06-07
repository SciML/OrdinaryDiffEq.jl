using ModelingToolkit, OrdinaryDiffEq, Test

@parameters t g e b
@variables v(t) w(t) F(t)
@derivatives D' ~ t
single_neuron_eqs = [
    D(v) ~ min(max(-2 - v, v), 2 - v) - w + F, # add the flux term
    D(w) ~ e * (v - g * w + b),
]
n1 = ODESystem(single_neuron_eqs, t, [v, w, F], [g, e, b], name = :n1)
n2 = ODESystem(single_neuron_eqs, t, [v, w, F], [g, e, b], name = :n2)
@parameters D Dk
connections = [0 ~ n1.F - D * Dk * max(n1.v - n2.v, 0)
    0 ~ n2.F - D * max(n2.v - n1.v, 0)]
connected = ODESystem(connections, t, [], [D, Dk], systems = [n1, n2], name = :connected)

u0 = [
    n1.v => -2,
    n1.w => -2 / 3,
    n1.F => 0,
    n2.v => -2,
    n2.w => -0.7,
    n2.F => 0,
]
tspan = (0.0, 1750.0)
p0 = [
    n1.g => 0.8,
    n1.e => 0.04,
    n1.b => 0,
    n2.g => 0.8,
    n2.e => 0.04,
    n2.b => 0.2,
    D => 0.047,
    Dk => 1,
]

prob = ODEProblem(connected, u0, tspan, p0)
sol = solve(prob, Rodas5(), initializealg = BrownFullBasicInit())
@test prob.u0 == sol[1]
sol = solve(prob, Rodas5(), initializealg = ShampineCollocationInit())
@test prob.u0 == sol[1]
