using StochasticDiffEq, Test, Random
using SDEProblemLibrary: prob_sde_linear
Random.seed!(100)
prob = prob_sde_linear

integrator = init(prob, EM(), dt = 1 // 2^(4), tstops = [0.33])

for (u, t) in tuples(integrator)
    @show u, t
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
