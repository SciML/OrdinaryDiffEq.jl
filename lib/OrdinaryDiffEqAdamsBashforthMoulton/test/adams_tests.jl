using OrdinaryDiffEqAdamsBashforthMoulton, DiffEqDevTools, DiffEqBase, Test
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

@testset "ABM steady-step function evaluations" begin
    for (alg, startup) in ((ABM32(), 2), (ABM43(), 3), (ABM54(), 4))
        for inplace in (false, true)
            calls = Ref(0)
            prob = if inplace
                ODEProblem(
                    (du, u, p, t) -> (calls[] += 1; du[1] = -u[1]),
                    [1.0], (0.0, 1.0)
                )
            else
                ODEProblem((u, p, t) -> (calls[] += 1; -u), 1.0, (0.0, 1.0))
            end
            integrator = init(prob, alg; dt = 1 / 32, adaptive = false)
            for _ in 1:startup
                step!(integrator)
            end
            for _ in 1:4
                nf = integrator.stats.nf
                old_calls = calls[]
                step!(integrator)
                @test integrator.stats.nf - nf == 2
                @test calls[] - old_calls == 2
            end
        end
    end
end

@testset "Fixed-step AB/ABM in-place vs out-of-place nf counts" begin
    function lorenz_iip!(du, u, p, t)
        p[] += 1
        du[1] = 10.0 * (u[2] - u[1])
        du[2] = u[1] * (28.0 - u[3]) - u[2]
        du[3] = u[1] * u[2] - (8 / 3) * u[3]
        return nothing
    end
    function lorenz_oop(u, p, t)
        p[] += 1
        return [
            10.0 * (u[2] - u[1]),
            u[1] * (28.0 - u[3]) - u[2],
            u[1] * u[2] - (8 / 3) * u[3],
        ]
    end
    u0 = [1.0, 0.0, 0.0]
    tspan = (0.0, 1.0)
    dt = 1 / 256
    for alg in (AB3(), AB4(), AB5(), ABM32(), ABM43(), ABM54())
        calls_iip = Ref(0)
        calls_oop = Ref(0)
        sol_iip = solve(
            ODEProblem(lorenz_iip!, copy(u0), tspan, calls_iip), alg; dt, adaptive = false
        )
        sol_oop = solve(
            ODEProblem(lorenz_oop, copy(u0), tspan, calls_oop), alg; dt, adaptive = false
        )
        @test sol_iip.stats.nf == sol_oop.stats.nf
        @test sol_iip.stats.nf == calls_iip[]
        @test sol_oop.stats.nf == calls_oop[]
    end
end

for i in 1:2
    prob = probArr[i]
    dt = 1 // 256
    integrator = init(prob, VCAB3(); dt, adaptive = false)
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

@testset "Float32 in-place Ralston start" begin
    f!(du, u, p, t) = (du .= -u; nothing)
    prob = ODEProblem(f!, Float32[1, 0.5], (0.0f0, 1.0f0))
    for alg in (AB3(), ABM32())
        sol = solve(prob, alg; dt = 0.01f0)
        @test eltype(sol.u[end]) == Float32
        @test sol.u[end] ≈ Float32[1, 0.5] * exp(-1.0f0) rtol = 1.0e-4
    end
end
