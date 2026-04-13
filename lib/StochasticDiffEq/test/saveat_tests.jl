using StochasticDiffEq, Test, Random
using SDEProblemLibrary: prob_sde_linear, prob_sde_2Dlinear
Random.seed!(100)
prob = prob_sde_linear
f(u, p, t) = 2u
prob = SDEProblem{false}(f, prob.f, prob.u0, prob.tspan)
sol = solve(prob, EM(), dt = 1 // 2^(4), saveat = 0.33)

sol.t == [0.0, 0.33, 0.33 + 0.33, 1]

sol = solve(prob, EM(), dt = 1 // 2^(4), saveat = 0.33, save_start = false)

sol.t == [0.33, 0.33 + 0.33, 1]
sol(0.4)

sol = solve(
    prob, EM(), dt = 1 // 2^(4), saveat = 0.33, save_start = false, save_everystep = true
)

@test 0.33 ∈ sol.t
@test 0.66 ∈ sol.t
sol.t != [0.33, 0.33 + 0.33, 1]

prob = prob_sde_2Dlinear
f(du, u, p, t) = prob.f(du, u, p, t)
prob2 = SDEProblem(f, prob.g, vec(prob.u0), prob.tspan)
Random.seed!(200)
sol = solve(
    prob2, EM(), dt = 1 // 2^(4), saveat = 0.33, save_start = false, save_everystep = true
)
Random.seed!(200)
sol2 = solve(
    prob2, EM(), dt = 1 // 2^(4), saveat = 0.33,
    save_start = false, save_everystep = true, save_idxs = 1:2:5
)

@test sol(0.43)[1:2:5] == sol2(0.43)

# Test save_on switch
sol = solve(
    prob2, EM(), dt = 1 // 2^(4), save_on = false, save_start = false, save_end = false
)
@test isempty(sol.t) && isempty(sol.u)
sol = solve(
    prob2, EM(), dt = 1 // 2^(4), saveat = 0.33,
    save_on = false, save_start = false, save_end = false
)
@test isempty(sol.t) && isempty(sol.u)

_saveat = [0.0, 0.25, 0.5, 1.0]
integ = init(SDEProblem((u, p, t) -> u, (u, p, t) -> u, 0.0, (0.0, 1.0)), SOSRI(), saveat = _saveat)
add_tstop!(integ, 2.0)
solve!(integ)
@test integ.sol.t == _saveat

integ = init(
    SDEProblem((u, p, t) -> u, (u, p, t) -> u, 0.0, (0.0, 1.0)), SOSRI(), saveat = _saveat, save_end = true
)
add_tstop!(integ, 2.0)
solve!(integ)
@test integ.sol.t == [0.0, 0.25, 0.5, 1.0, 2.0]

integ = init(
    SDEProblem((u, p, t) -> u, (u, p, t) -> u, 0.0, (0.0, 1.0)), SOSRI(), saveat = _saveat, save_end = false
)
add_tstop!(integ, 2.0)
solve!(integ)
@test integ.sol.t == _saveat

# Catch save for maxiters
ode = SDEProblem((u, p, t) -> u, (u, p, t) -> u, 1.0, (0.0, 1.0))
sol = solve(ode, SOSRI(), save_everystep = false) # okay, as expected
@test length(sol) == 2
@info "Warning Expected"
sol = solve(ode, SOSRI(), save_everystep = false, maxiters = 3) # doesn't save the final solution anymore!
@test length(sol) == 2
