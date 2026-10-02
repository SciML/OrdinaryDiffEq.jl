using OrdinaryDiffEqBDF, LinearAlgebra, Test

const A = [-1.0 0.35 0.1; 0.2 -1.4 0.25; 0.15 0.3 -1.2]
const U0 = [1.0, -0.5]
const UNEW = 0.3
const TG = 0.5
const TS = 1.0
const TF = 2.0

f!(du, u, p, t) = (mul!(du, view(A, 1:length(u), 1:length(u)), u); nothing)

function exact(t)
    u = exp(A[1:2, 1:2] * min(t, TG)) * U0
    t <= TG && return u
    u3 = exp(A * (min(t, TS) - TG)) * vcat(u, UNEW)
    t <= TS && return u3
    return exp(A[1:2, 1:2] * (t - TS)) * u3[1:2]
end

function resize_callback!(disc, restart_orders)
    grew = Ref(false)
    shrank = Ref(false)
    function condition(u, t, integrator)
        return (!grew[] && t >= TG) || (!shrank[] && t >= TS)
    end
    function affect(integrator)
        n = grew[] ? 2 : 3
        resize!(integrator, n)
        n == 3 && (integrator.u[3] = UNEW)
        push!(restart_orders, integrator.cache.order)
        grew[] ? (shrank[] = true) : (grew[] = true)
        return derivative_discontinuity!(integrator, disc)
    end
    return DiscreteCallback(condition, affect; save_positions = (false, false))
end

@testset "BDF state resizing" begin
    prob = ODEProblem(f!, copy(U0), (0.0, TF))

    for alg in (FBDF(), QNDF()), disc in (true, false)
        restart_orders = Int[]
        sol = solve(
            prob, alg; callback = resize_callback!(disc, restart_orders),
            tstops = [TG, TS], reltol = 1.0e-10, abstol = 1.0e-12
        )
        @test sol(0.75) ≈ exact(0.75) rtol = 1.0e-7 atol = 1.0e-10
        @test sol.u[end] ≈ exact(TF) rtol = 1.0e-7 atol = 1.0e-10
        @test restart_orders == [1, 1]
    end
end
