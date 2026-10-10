using DelayDiffEq, OrdinaryDiffEqRosenbrock
using SciMLBase: ReturnCode, first_tstop
using Test

# A discontinuity found by tracking a dependent delay lies one ulp after another
# discontinuity. Handling the first one shifts `t` by one ulp onto the tracked one,
# which must be dropped with it rather than left as a stop at the current time,
# where it would force a zero-length step (and a non-finite Rosenbrock W).
@testset "tracked discontinuity one ulp after a handled discontinuity (neutral = $neutral)" for neutral in (false, true)
    f!(du, u, h, p, t) = (du[1] = -h(p, t - 1)[1]; nothing)
    prob = DDEProblem(
        f!, [1.0], (p, t) -> [1.0], (0.0, 2.0);
        dependent_lags = ((u, p, t) -> 1.0,), neutral
    )
    integrator = init(prob, MethodOfSteps(Rodas5P()); dt = 0.1)
    step!(integrator)

    # track the propagated discontinuity as a rejected step over [t, 1.5] would
    dt = integrator.dt
    integrator.dt = 1.5 - integrator.t
    DelayDiffEq.track_propagated_discontinuities!(integrator)
    integrator.dt = dt
    td = first_tstop(integrator)
    @test td ≈ 1

    push!(integrator.opts.d_discontinuities, Discontinuity(prevfloat(td), 1))
    add_tstop!(integrator, prevfloat(td))

    sol = solve!(integrator)
    @test sol.retcode == ReturnCode.Success
    @test prevfloat(td) in sol.t
    @test sol.t[end] == 2
end
