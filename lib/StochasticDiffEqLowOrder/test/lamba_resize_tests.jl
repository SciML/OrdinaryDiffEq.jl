using StochasticDiffEqLowOrder
using Test

function f_diag!(du, u, p, t)
    @. du = 1.01 * u
    return nothing
end

function g_diag!(du, u, p, t)
    @. du = 0.87 * u
    return nothing
end

@testset "LambaEM/LambaEulerHeun DiscreteCallback resize! (diagonal noise)" begin
    for alg in (LambaEM(), LambaEulerHeun())
        u0 = ones(4)
        tspan = (0.0, 1.0)
        prob = SDEProblem(f_diag!, g_diag!, u0, tspan)
        fired = Ref(false)
        condition = (u, t, integrator) -> !fired[] && t >= 0.5
        affect! = function (integrator)
            fired[] = true
            n = length(integrator.u)
            resize!(integrator, n + 2)
            integrator.u[(n + 1):end] .= 1.0
            return nothing
        end
        cb = DiscreteCallback(condition, affect!)
        sol = solve(
            prob, alg; dt = 0.01, adaptive = false, dense = false,
            callback = cb, save_everystep = false
        )
        @test fired[]
        @test length(sol.u[end]) == 6
        @test all(isfinite, sol.u[end])
    end
end
