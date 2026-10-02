using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, Test

f!(du, u, p, t) = (du .= -u)

function resize_callback!(restart_orders)
    grew = Ref(false)
    shrank = Ref(false)
    function condition(u, t, integrator)
        return (!grew[] && t >= 0.5) || (!shrank[] && t >= 1.0)
    end
    function affect(integrator)
        if !grew[]
            resize!(integrator, 3)
            integrator.u[3] = 0.5
            grew[] = true
        else
            resize!(integrator, 2)
            shrank[] = true
        end
        hasproperty(integrator.cache, :order) && push!(restart_orders, integrator.cache.order)
        return derivative_discontinuity!(integrator, true)
    end
    return DiscreteCallback(condition, affect; save_positions = (false, false))
end

@testset "BDF state resizing" begin
    prob = ODEProblem(f!, ones(2), (0.0, 2.0))
    reference = solve(
        prob, ImplicitEuler(); callback = resize_callback!(Int[]),
        tstops = [0.5, 1.0], reltol = 1.0e-9, abstol = 1.0e-11
    )

    for alg in (FBDF(), QNDF())
        restart_orders = Int[]
        sol = solve(
            prob, alg; callback = resize_callback!(restart_orders),
            tstops = [0.5, 1.0], reltol = 1.0e-9, abstol = 1.0e-11
        )
        @test length(sol.u[end]) == 2
        @test sol.u[end] ≈ reference.u[end] rtol = 5.0e-5 atol = 2.0e-7
        @test restart_orders == [1, 1]
    end
end
