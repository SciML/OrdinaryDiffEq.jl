using StochasticDiffEq, DiffEqNoiseProcess, Test

println("Scalar g")
A = [-1.0 0.0; 0.0 -0.5]
u0 = [1.0, 1.0];
tspan = (0.0, 1.0)
_f = (u, p, t) -> t * (A * u)
_g = (u, p, t) -> 1.0
prob = SDEProblem(SDEFunction(_f, _g), u0, tspan)
integrator = init(prob, SKenCarp(); adaptive = false, dt = 0.01)
step!(integrator)
@test_broken solve(prob, SOSRI(); adaptive = false, dt = 0.01) isa RODESolution
solve(prob, SOSRA(); adaptive = false, dt = 0.01)
solve(prob, EM(); adaptive = false, dt = 0.01)
@test_broken solve(prob, RKMil(); adaptive = false, dt = 0.01) isa RODESolution
solve(prob, SRIW1(); adaptive = false, dt = 0.01)
@test_broken solve(prob, SOSRI2(); adaptive = false, dt = 0.01) isa RODESolution
solve(prob, SOSRA2(); adaptive = false, dt = 0.01)

println("Vector g")
_g = (u, p, t) -> [1.0, 1.0]
prob = SDEProblem(SDEFunction(_f, _g), u0, tspan)
println("Implicit EM")
integrator = init(prob, ImplicitEM(); adaptive = false, dt = 0.01)
step!(integrator)
println("SKenCarp")
integrator = init(prob, SKenCarp(); adaptive = false, dt = 0.01)
step!(integrator)
solve(prob, SOSRI(); adaptive = false, dt = 0.01)
solve(prob, SOSRA(); adaptive = false, dt = 0.01)
solve(prob, EM(); adaptive = false, dt = 0.01)
solve(prob, RKMil(); adaptive = false, dt = 0.01)
solve(prob, SOSRI2(); adaptive = false, dt = 0.01)
solve(prob, SOSRA2(); adaptive = false, dt = 0.01)
