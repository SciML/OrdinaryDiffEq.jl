using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase, Test
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear

probArr = Vector{ODEProblem}(undef, 2)

probArr[1] = prob_ode_linear
probArr[2] = prob_ode_2Dlinear

function fixed_step_ϕstar(k)
    ∇ = Vector{typeof(k[end][1])}(undef, 3)
    ∇[1] = k[end][1]
    ∇[2] = ∇[1] - k[end - 1][1]
    ∇[3] = ∇[2] - k[end - 1][1] + k[end - 2][1]
    return ∇
end

for i in 1:2
    prob = probArr[i]
    dt = 1 // 256
    integrator = init(prob, VCAB3(), dt = dt, adaptive = false)
    for i in 1:3
        step!(integrator)
    end
    @test integrator.cache.g == [1, 1 / 2, 5 / 12] * dt
    step!(integrator)
    @test integrator.cache.g == [1, 1 / 2, 5 / 12] * dt
    step!(integrator)
    @test integrator.cache.g == [1, 1 / 2, 5 / 12] * dt
end

for i in 1:2
    prob = probArr[i]
    # VCAB3
    integrator = init(prob, VCAB3(), dt = 1 // 256, adaptive = false)
    for i in 1:3
        step!(integrator)
    end
    # in perform_step, after swapping array using pointer, ϕstar_nm1 points to ϕstar_n
    @test integrator.cache.ϕstar_nm1 == fixed_step_ϕstar(integrator.sol.k)
    step!(integrator)
    @test integrator.cache.ϕstar_nm1 == fixed_step_ϕstar(integrator.sol.k)
    step!(integrator)
    @test integrator.cache.ϕstar_nm1 == fixed_step_ϕstar(integrator.sol.k)

    # VCAB4
    sol1 = solve(prob, VCAB4(), dt = 1 // 256, adaptive = false)
    sol2 = solve(prob, AB4(), dt = 1 // 256)
    @test sol1.u ≈ sol2.u

    # VCAB5
    sol1 = solve(prob, VCAB5(), dt = 1 // 256, adaptive = false)
    sol2 = solve(prob, AB5(), dt = 1 // 256)
    @test sol1.u ≈ sol2.u
end