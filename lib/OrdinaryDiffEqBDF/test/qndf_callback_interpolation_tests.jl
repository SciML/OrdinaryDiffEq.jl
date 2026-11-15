using OrdinaryDiffEqBDF, Test

# A ContinuousCallback truncates the step it finds its event in, and the step saved up to
# the event keeps QNDF's backward differences for interpolation. They must describe the
# completed step's polynomial on the truncated step (rebased to its length and endpoint):
# the dense output inside the truncated step is then as accurate as elsewhere.

# u = (sin t, cos t, exp(-t)): smooth, with the event where sin t crosses 0.9
f!(du, u, p, t) = (du[1] = u[2]; du[2] = -u[1]; du[3] = -u[3]; nothing)
f(u, p, t) = [u[2], -u[1], -u[3]]
exact(t) = [sin(t), cos(t), exp(-t)]
condition(u, t, integrator) = u[1] - 0.9
affect!(integrator) = nothing

# The largest error of the dense output inside the steps that end at the event, and
# inside all other steps
function dense_errors(prob, alg, t_event)
    sol = solve(prob, alg; callback = ContinuousCallback(condition, affect!), abstol = 1.0e-10, reltol = 1.0e-10)
    at_event, elsewhere = 0.0, 0.0
    for k in 1:(length(sol.t) - 1)
        a, b = sol.t[k], sol.t[k + 1]
        b > a || continue
        e = maximum(maximum(abs.(sol(τ) .- exact(τ))) for τ in range(a, b, length = 21))
        if abs(b - t_event) < 1.0e-8
            at_event = max(at_event, e)
        else
            elsewhere = max(elsewhere, e)
        end
    end
    return at_event, elsewhere
end

t_event = asin(0.9)
@testset "QNDF dense output in a step truncated by a callback" begin
    for prob in (ODEProblem(f!, [0.0, 1.0, 1.0], (0.0, 2.0)), ODEProblem(f, [0.0, 1.0, 1.0], (0.0, 2.0)))
        for alg in (QNDF(), QBDF(), FBDF())
            at_event, elsewhere = dense_errors(prob, alg, t_event)
            @test at_event > 0              # the event's step was found
            @test at_event < 10 * max(elsewhere, 1.0e-9)
        end
    end
end

# The scalar path (the constant cache with a Number state): u = exp(-t), event at 0.5
@testset "QNDF dense output in a truncated step, scalar state" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 2.0))
    t_half = log(2.0)
    sol = solve(
        prob, QNDF(); callback = ContinuousCallback((u, t, integrator) -> u - 0.5, affect!),
        abstol = 1.0e-10, reltol = 1.0e-10
    )
    k = findfirst(t -> abs(t - t_half) < 1.0e-8, sol.t)
    @test k !== nothing
    a, b = sol.t[k - 1], sol.t[k]
    @test maximum(abs(sol(τ) - exp(-τ)) for τ in range(a, b, length = 21)) < 1.0e-8
end
