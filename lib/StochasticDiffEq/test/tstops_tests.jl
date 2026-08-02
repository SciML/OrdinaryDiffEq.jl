using StochasticDiffEq, Test, Random
using DiffEqNoiseProcess: WienerProcess, RSWM
using SDEProblemLibrary: prob_sde_linear
Random.seed!(100)
prob = prob_sde_linear

integrator = init(prob, EM(), dt = 1 // 2^(4), tstops = [0.33])

# SciMLBase v3: `tuples(integrator)` removed; iterate integrator directly.
for integ in integrator
    @show integ.u, integ.t
end

sol = solve(prob, EM(), dt = 1 // 2^(4), tstops = [0.33])

@test 0.33 ∈ sol.t

sol = solve(prob, EM(), tstops = [0.33, 0.8, 1.0])

@test sol.t == [0.0, 0.33, 0.8, 1.0]

sol = solve(prob, SRIW1(), tstops = [0.33])

@test 0.33 ∈ sol.t

# check reverse time and negative start times
for (i, tdir) in enumerate([-1.0; 1.0])
    @info i
    prob2 = remake(prob_sde_linear, tspan = (tdir * 1.0, 0.0))
    integrator = init(prob2, SRIW1())
    tstops = tdir .* [0, 0.33, 0.8, 1]
    for tstop in tstops
        add_tstop!(integrator, tstop)
    end
    solve!(integrator)
    for tstop in tstops
        @test tstop ∈ integrator.sol.t
    end
end

# SciML/OrdinaryDiffEq.jl#3175 / StochasticDiffEq.jl#413:
# late add_tstop! after init must shorten the first EM step to the tstop.
f3175(u, p, t) = 1.0
g3175(u, p, t) = 0.1
prob3175 = SDEProblem(f3175, g3175, [1.0], (0.0, 1.0))

i_em = init(prob3175, EM(); dt = 0.02)
add_tstop!(i_em, 0.01)
step!(i_em)
@test i_em.t == 0.01

i_lamba = init(prob3175, LambaEM(); dt = 0.02)
add_tstop!(i_lamba, 0.01)
step!(i_lamba)
@test i_lamba.t == 0.01

# Past-time add_tstop! is rejected for both fixed-step and adaptive EM.
i_past = init(prob3175, EM(); dt = 0.02)
step!(i_past)
@test_throws ErrorException add_tstop!(i_past, 0.01)

# The noise has to follow the solver onto the shortened step. Without bridging,
# W keeps stepping on the original dt grid and drifts away from sol.t for the
# rest of the integration.
for alg in (EM(), LambaEM())
    i = init(prob3175, alg; dt = 0.02, save_noise = true)
    add_tstop!(i, 0.01)
    solve!(i)
    @test i.W.t == i.sol.t
    @test i.W.curt == i.t
end

# ... and bridging is what it must do: committing the pending increment instead
# would make the noise take a step the solver never took.
for alg in (EM(), LambaEM())
    i = init(prob3175, alg; dt = 0.02)
    add_tstop!(i, 0.01)
    solve!(i)
    @test i.W.iter == i.iter
end

# For du = dW the Euler step is exact, so the endpoints must agree with the
# Brownian path the solution was built from.
fexact(u, p, t) = 0.0
gexact(u, p, t) = 1.0
probexact = SDEProblem(fexact, gexact, [0.0], (0.0, 1.0))

for alg in (EM(), LambaEM())
    i = init(probexact, alg; dt = 0.02, save_noise = true)
    add_tstop!(i, 0.01)
    solve!(i)
    @test i.sol.u[end][1] ≈ only(i.W.curW) atol = 1.0e-12
end

# RSwM1 hands a stack chunk back through W.dt, so the tstop shortening has to
# bridge that chunk rather than relabel it.
for adaptivealg in (:RSwM1, :RSwM2, :RSwM3)
    W1 = WienerProcess(0.0, 0.0; rswm = RSWM(adaptivealg = adaptivealg))
    prob1 = SDEProblem(fexact, gexact, 0.0, (0.0, 1.0), noise = W1)
    i = init(prob1, LambaEM(); dt = 0.02, tstops = [0.331], save_noise = true)
    add_tstop!(i, 0.01)
    solve!(i)
    @test 0.01 ∈ i.sol.t
    @test 0.331 ∈ i.sol.t
    @test i.sol.u[end] ≈ i.W.curW atol = 1.0e-12
end
