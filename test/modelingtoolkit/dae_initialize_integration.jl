using ModelingToolkit, OrdinaryDiffEq, NonlinearSolve, Test
using ModelingToolkit: D_nounits as D, t_nounits as t

@parameters g e b
@variables v(t) w(t) F(t)
single_neuron_eqs = [
    D(v) ~ min(max(-2 - v, v), 2 - v) - w + F, # add the flux term
    D(w) ~ e * (v - g * w + b),
]
n1 = System(single_neuron_eqs, t, [v, w, F], [g, e, b], name = :n1)
n2 = System(single_neuron_eqs, t, [v, w, F], [g, e, b], name = :n2)
@parameters Di Dk
connections = [
    0 ~ n1.F - Di * Dk * max(n1.v - n2.v, 0)
    0 ~ n2.F - Di * max(n2.v - n1.v, 0)
]
connected = System(connections, t, [], [Di, Dk], systems = [n1, n2], name = :connected)
connected = complete(connected)

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
    Di => 0.047,
    Dk => 1,
]

prob = ODEProblem(connected, [u0; p0], tspan)
sol = solve(prob, Rodas5(), initializealg = BrownFullBasicInit())
@test prob.u0 == sol.u[1]
sol = solve(prob, Rodas5(), initializealg = ShampineCollocationInit())
@test prob.u0 == sol.u[1]
#test initialization when given a specific nonlinear solver
using NonlinearSolve
sol = solve(prob, Rodas5(), initializealg = BrownFullBasicInit(1.0e-10, RobustMultiNewton()))
@test prob.u0 == sol.u[1]

# Initialize on ODEs
# https://github.com/SciML/ModelingToolkit.jl/issues/2508

function testsys(du, u, p, t)
    return du[1] = -2
end
function initsys(du, u, p)
    return du[1] = -1 + u[1]
end
nlprob = NonlinearProblem(initsys, [0.0])
initprobmap(nlprob) = nlprob.u
sol = solve(nlprob)

_f = ODEFunction(testsys; initializeprob = nlprob, initializeprobmap = initprobmap)
prob = ODEProblem(_f, [0.0], (0.0, 1.0))
sol = solve(prob, Tsit5())
@test SciMLBase.successful_retcode(sol)
@test sol.u[1] == [1.0]

prob = ODEProblem(_f, [0.0], (0.0, 1.0))
sol = solve(prob, Tsit5(), dt = 1.0e-10)
@test SciMLBase.successful_retcode(sol)
@test sol.u[1] == [1.0]
@test sol.u[2] ≈ [0.9999999998]
@test sol.u[end] ≈ [-1.0]

sol = solve(prob, Rodas5P(), dt = 1.0e-10)
@test SciMLBase.successful_retcode(sol)
@test sol.u[1] == [1.0]
@test sol.u[2] ≈ [0.9999999998]
@test sol.u[end] ≈ [-1.0]

@testset "`reinit!` updates initial parameters" begin
    # https://github.com/SciML/ModelingToolkit.jl/issues/3451
    # https://github.com/SciML/ModelingToolkit.jl/issues/3504
    @variables x(t) y(t)
    @parameters c1 c2
    @mtkcompile sys = System([D(x) ~ -c1 * x + c2 * y, D(y) ~ c1 * x - c2 * y], t)
    prob = ODEProblem(sys, [x => 1.0, y => 2.0, c1 => 1.0, c2 => 2.0], (0.0, 1.0))
    @test prob.ps[Initial(x)] ≈ 1.0
    @test prob.ps[Initial(y)] ≈ 2.0
    integ = init(prob, Tsit5())
    @test integ.ps[Initial(x)] ≈ 1.0
    @test integ.ps[Initial(y)] ≈ 2.0
    new_u0 = ModelingToolkit.get_u0(sys, Dict(x => 2.0, y => 3.0))
    reinit!(integ, new_u0)
    @test integ.ps[Initial(x)] ≈ 2.0
    @test integ.ps[Initial(y)] ≈ 3.0
    @test integ[x] ≈ 2.0
    @test integ[y] ≈ 3.0
end
