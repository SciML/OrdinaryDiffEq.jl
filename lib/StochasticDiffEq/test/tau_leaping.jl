using StochasticDiffEq, JumpProcesses, DiffEqBase, Statistics, OrdinaryDiffEq
using Test, LinearAlgebra

function regular_rate(out, u, p, t)
    out[1] = (0.1 / 1000.0) * u[1] * u[2]
    return out[2] = 0.01u[2]
end

const dc = zeros(3, 2)
dc[1, 1] = -1
dc[2, 1] = 1
dc[2, 2] = -1
dc[3, 2] = 1

function regular_c(du, u, p, t, counts, mark)
    return mul!(du, dc, counts)
end

rj = RegularJump(regular_rate, regular_c, 2)
jumps = JumpSet(rj)
iip_prob = DiscreteProblem([999.0, 1, 0], (0.0, 250.0))
jump_iipprob = JumpProblem(iip_prob, Direct(), rj)
jump_iipprob_pure = JumpProblem(iip_prob, PureLeaping(), rj)
@time sol = solve(jump_iipprob, TauLeaping())
@time sol = solve(jump_iipprob_pure, SimpleTauLeaping(); dt = 1.0)
@time sol = solve(jump_iipprob, TauLeaping(); dt = 1.0, adaptive = false)
@time sol = solve(jump_iipprob, CaoTauLeaping(); dt = 1.0)
@time sol = solve(jump_iipprob, CaoTauLeaping())

N = 40_000
sol1 = solve(EnsembleProblem(jump_iipprob_pure), SimpleTauLeaping(); dt = 1.0, trajectories = N)
sol2 = solve(
    EnsembleProblem(jump_iipprob), TauLeaping(); dt = 1.0,
    adaptive = false, save_everystep = false, trajectories = N
)

mean1 = mean([sol1.u[i][end, end] for i in 1:N])
mean2 = mean([sol2.u[i][end, end] for i in 1:N])
@test mean1 ≈ mean2 rtol = 1.0e-2

f(du, u, p, t) = (du .= 0)
g(du, u, p, t) = (du .= 0)
iip_sdeprob = SDEProblem(f, g, [999.0, 1, 0], (0.0, 250.0))
jumpdiff_iipprob = JumpProblem(iip_sdeprob, Direct(), rj)
@time sol = solve(jumpdiff_iipprob, EM(); dt = 1.0)
@time sol = solve(jumpdiff_iipprob, ImplicitEM(); dt = 1.0, adaptive = false)

sol = solve(EnsembleProblem(jumpdiff_iipprob), EM(); dt = 1.0, trajectories = N)
meanX = mean([sol.u[i][end, end] for i in 1:N])
@test mean1 ≈ meanX rtol = 1.0e-2

sol = solve(EnsembleProblem(jumpdiff_iipprob), ImplicitEM(); dt = 1.0, trajectories = N)
meanX = mean([sol.u[i][end, end] for i in 1:N])
@test mean1 ≈ meanX rtol = 1.0e-2

iip_prob = DiscreteProblem([999, 1, 0], (0.0, 250.0))
jump_iipprob = JumpProblem(iip_prob, Direct(), rj)
sol = solve(jump_iipprob, TauLeaping())

function rate_oop(u, p, t)
    return [(0.1 / 1000.0) * u[1] * u[2], 0.01u[2]]
end

function regular_c(u, p, t, counts, mark)
    return dc * counts
end

rj = RegularJump(rate_oop, regular_c, 2)
jumps = JumpSet(rj)
prob = DiscreteProblem([999.0, 1, 0], (0.0, 250.0))
jump_prob = JumpProblem(prob, Direct(), rj)
sol = solve(jump_prob, TauLeaping(), reltol = 5.0e-2)

sol2 = solve(
    EnsembleProblem(jump_prob), TauLeaping(); dt = 1.0,
    adaptive = false, save_everystep = false, trajectories = N
)
mean2 = mean([sol2.u[i][end, end] for i in 1:N])
@test mean1 ≈ mean2 rtol = 1.0e-2

foop(u, p, t) = [0.0, 0.0, 0.0]
goop(u, p, t) = [0.0, 0.0, 0.0]
oop_sdeprob = SDEProblem(foop, goop, [999.0, 1, 0], (0.0, 250.0))
jumpdiff_prob = JumpProblem(oop_sdeprob, Direct(), rj)
@time sol = solve(jumpdiff_prob, EM(); dt = 1.0)
@time sol = solve(jumpdiff_prob, ImplicitEM(); dt = 1.0)

sol = solve(EnsembleProblem(jumpdiff_prob), EM(); dt = 1.0, trajectories = 10_000)
meanX = mean([sol.u[i][end, end] for i in 1:10_000])
@test mean1 ≈ meanX rtol = 1.0e-2

sol = solve(EnsembleProblem(jumpdiff_prob), ImplicitEM(); dt = 1.0, trajectories = 1_000)
meanX = mean([sol.u[i][end, end] for i in 1:1_000])
@test mean1 ≈ meanX rtol = 1.0e-1

# ThetaTrapezoidalTauLeaping: Implicit weak second order tau-leaping method
# Test in-place version with fixed dt
function regular_rate_iip(out, u, p, t)
    out[1] = (0.1 / 1000.0) * u[1] * u[2]
    return out[2] = 0.01u[2]
end

const dc_iip = zeros(3, 2)
dc_iip[1, 1] = -1
dc_iip[2, 1] = 1
dc_iip[2, 2] = -1
dc_iip[3, 2] = 1

function regular_c_iip(du, u, p, t, counts, mark)
    return mul!(du, dc_iip, counts)
end

rj_iip = RegularJump(regular_rate_iip, regular_c_iip, 2)
iip_prob_theta = DiscreteProblem([999.0, 1, 0], (0.0, 250.0))
jump_iipprob_theta = JumpProblem(iip_prob_theta, Direct(), rj_iip)

# Test basic solve
@time sol_theta = solve(
    jump_iipprob_theta, ThetaTrapezoidalTauLeaping(); dt = 1.0, adaptive = false
)

# Test with different theta values
for theta in [0.25, 0.5, 0.75]
    @time sol = solve(
        jump_iipprob_theta, ThetaTrapezoidalTauLeaping(; theta = theta);
        dt = 1.0, adaptive = false
    )
    @test length(sol.t) > 0
end

# Compare mean with TauLeaping - should give similar results
N_theta = 10_000
sol_tauleaping_theta = solve(
    EnsembleProblem(jump_iipprob_theta), TauLeaping(); dt = 1.0,
    adaptive = false, save_everystep = false, trajectories = N_theta
)
sol_theta_ens = solve(
    EnsembleProblem(jump_iipprob_theta), ThetaTrapezoidalTauLeaping();
    dt = 1.0, adaptive = false, save_everystep = false, trajectories = N_theta
)

mean_tauleaping_theta = mean([sol_tauleaping_theta.u[i][end, end] for i in 1:N_theta])
mean_theta_ens = mean([sol_theta_ens.u[i][end, end] for i in 1:N_theta])
@test mean_tauleaping_theta ≈ mean_theta_ens rtol = 5.0e-2

# Test out-of-place version
function rate_oop_theta(u, p, t)
    return [(0.1 / 1000.0) * u[1] * u[2], 0.01u[2]]
end

const dc_oop_theta = zeros(3, 2)
dc_oop_theta[1, 1] = -1
dc_oop_theta[2, 1] = 1
dc_oop_theta[2, 2] = -1
dc_oop_theta[3, 2] = 1

function regular_c_oop_theta(u, p, t, counts, mark)
    return dc_oop_theta * counts
end

rj_oop = RegularJump(rate_oop_theta, regular_c_oop_theta, 2)
prob_oop_theta = DiscreteProblem([999.0, 1, 0], (0.0, 250.0))
jump_prob_oop_theta = JumpProblem(prob_oop_theta, Direct(), rj_oop)

@time sol_oop = solve(
    jump_prob_oop_theta, ThetaTrapezoidalTauLeaping(); dt = 1.0, adaptive = false
)
@test length(sol_oop.t) > 0

# Ensemble test for out-of-place
sol_theta_oop_ens = solve(
    EnsembleProblem(jump_prob_oop_theta),
    ThetaTrapezoidalTauLeaping();
    dt = 1.0, adaptive = false, save_everystep = false,
    trajectories = N_theta
)
mean_theta_oop = mean([sol_theta_oop_ens.u[i][end, end] for i in 1:N_theta])
@test mean_theta_oop ≈ mean_theta_ens rtol = 5.0e-2

# ImplicitTauLeaping: First-order implicit (backward Euler) tau-leaping method
# Reuse the same problem setup from ThetaTrapezoidalTauLeaping tests

# Test basic solve with in-place functions
@time sol_implicit = solve(
    jump_iipprob_theta, ImplicitTauLeaping(); dt = 1.0, adaptive = false
)
@test length(sol_implicit.t) > 0

# Test with out-of-place functions
@time sol_implicit_oop = solve(
    jump_prob_oop_theta, ImplicitTauLeaping(); dt = 1.0, adaptive = false
)
@test length(sol_implicit_oop.t) > 0

# Compare mean with TauLeaping - should give similar results
N_implicit = 10_000
sol_implicit_ens = solve(
    EnsembleProblem(jump_iipprob_theta), ImplicitTauLeaping();
    dt = 1.0, adaptive = false, save_everystep = false, trajectories = N_implicit
)
mean_implicit_ens = mean([sol_implicit_ens.u[i][end, end] for i in 1:N_implicit])
@test mean_tauleaping_theta ≈ mean_implicit_ens rtol = 5.0e-2
