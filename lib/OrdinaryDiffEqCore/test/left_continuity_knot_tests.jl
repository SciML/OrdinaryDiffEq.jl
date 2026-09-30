using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, Test
using SciMLBase: u_modified!

# Regression for duplicated knots from save_positions / initialize callbacks:
# sol(t; continuity=:left) must return the pre-jump state even when saveat
# disables dense output (zero-width interval Θ fallback).

function make_jump_callback(t_event)
    affect!(integrator) = (integrator.u[1] += 1.0)
    return DiscreteCallback(
        (u, t, integrator) -> t == t_event,
        affect!;
        initialize = (c, u, t, integrator) -> (affect!(integrator); u_modified!(integrator, true)),
        save_positions = (true, true),
    )
end

@testset "left continuity at duplicated saveat knots (tdir>0)" begin
    cb = make_jump_callback(1.0)
    prob = ODEProblem((u, p, t) -> -u, [0.0], (0.0, 2.0))
    sol = solve(prob, Tsit5(); callback = cb, tstops = [1.0], saveat = [0.0, 1.0, 2.0])

    @test sol.t[1] == sol.t[2] == 0.0
    @test sol.u[1] == [0.0]
    @test sol.u[2] == [1.0]

    @test sol(0.0; continuity = :left) == [0.0]
    @test sol(0.0; continuity = :right) == [1.0]
    @test sol([0.0]; continuity = :left).u[1] == [0.0]
    @test sol([0.0]; continuity = :right).u[1] == [1.0]

    out = zeros(1)
    sol(out, 0.0; continuity = :left)
    @test out == [0.0]
    sol(out, 0.0; continuity = :right)
    @test out == [1.0]

    # Mid-interval callback knot still distinguishes pre/post
    i1 = findfirst(==(1.0), sol.t)
    @test sol.t[i1] == sol.t[i1 + 1] == 1.0
    @test sol(1.0; continuity = :left) == sol.u[i1]
    @test sol(1.0; continuity = :right) == sol.u[i1 + 1]
end

@testset "left continuity at duplicated knots with dense output" begin
    cb = make_jump_callback(1.0)
    prob = ODEProblem((u, p, t) -> -u, [0.0], (0.0, 2.0))
    sol = solve(prob, Tsit5(); callback = cb, tstops = [1.0])

    @test sol.t[1] == sol.t[2] == 0.0
    @test sol(0.0; continuity = :left) == [0.0]
    @test sol(0.0; continuity = :right) == [1.0]
end

@testset "left continuity at duplicated saveat knots (tdir<0)" begin
    cb = make_jump_callback(1.0)
    prob = ODEProblem((u, p, t) -> -u, [0.0], (2.0, 0.0))
    sol = solve(prob, Tsit5(); callback = cb, tstops = [1.0], saveat = [2.0, 1.0, 0.0])

    @test sol.t[1] == sol.t[2] == 2.0
    @test sol.u[1] == [0.0]
    @test sol.u[2] == [1.0]

    @test sol(2.0; continuity = :left) == [0.0]
    @test sol(2.0; continuity = :right) == [1.0]
    @test sol([2.0]; continuity = :left).u[1] == [0.0]

    out = zeros(1)
    sol(out, 2.0; continuity = :left)
    @test out == [0.0]
    sol(out, 2.0; continuity = :right)
    @test out == [1.0]
end
