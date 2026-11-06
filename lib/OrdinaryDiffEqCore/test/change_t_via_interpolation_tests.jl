using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, Test
using SciMLBase: change_t_via_interpolation!, successful_retcode

decay = ODEProblem((du, u, p, t) -> (du[1] = -u[1]; nothing), [1.0], (0.0, 4.0))

rewind_to_midstep(save_positions) = DiscreteCallback(
    (u, t, integrator) -> t == 2.0,
    function (integrator)
        change_t_via_interpolation!(
            integrator, integrator.tprev + (integrator.t - integrator.tprev) / 2
        )
        return nothing
    end;
    save_positions,
)

@testset "moving t behind a saved point is rejected" begin
    # `save_positions[1]` stores the step end before `affect!` runs, so the rewind would
    # put `integrator.t` behind `sol.t[end]`.
    @test_throws ErrorException solve(
        decay, Tsit5(); tstops = [2.0], callback = rewind_to_midstep((true, true))
    )
end

@testset "a rewind with nothing saved after it is left alone" begin
    sol = solve(
        decay, Tsit5(); tstops = [2.0], callback = rewind_to_midstep((false, false)),
        save_everystep = false
    )
    @test successful_retcode(sol)
    @test issorted(sol.t)
end

@testset "the root find of a continuous callback still rewinds" begin
    halved = ContinuousCallback((u, t, integrator) -> u[1] - 0.5, integrator -> nothing)
    sol = solve(decay, Tsit5(); callback = halved, abstol = 1.0e-10, reltol = 1.0e-10)
    @test successful_retcode(sol)
    @test issorted(sol.t)
    @test any(u -> isapprox(u[1], 0.5; atol = 1.0e-8), sol.u)
end
