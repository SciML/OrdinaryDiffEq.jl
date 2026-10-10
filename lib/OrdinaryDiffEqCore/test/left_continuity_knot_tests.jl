using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, Test
using SciMLBase: u_modified!

# Regression for duplicated knots from save_positions / initialize callbacks:
# sol(t; continuity=:left) must return the pre-jump state.

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

# Guard only: Tsit5 dense already returned the pre-jump state on master.
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

# Stiff dense interpolants divide by dt; zero-width intervals must collapse.
@testset "left continuity Rodas5P dense at duplicated initial knot" begin
    cb = make_jump_callback(1.0)
    prob = ODEProblem((u, p, t) -> -u, [0.0], (0.0, 2.0))
    sol = solve(prob, Rodas5P(); callback = cb, tstops = [1.0])

    @test sol.t[1] == sol.t[2] == 0.0
    @test sol.u[1] == [0.0]
    @test sol.u[2] == [1.0]

    left = sol(0.0; continuity = :left)
    @test left == [0.0]
    @test all(!isnan, left)
    @test sol(0.0; continuity = :right) == [1.0]

    left_vec = sol([0.0]; continuity = :left)
    @test left_vec.u[1] == [0.0]
    @test all(!isnan, left_vec.u[1])

    @test sol(0.0; idxs = 1, continuity = :left) == 0.0
    @test !isnan(sol(0.0; idxs = 1, continuity = :left))
end

# Extrapolation past a duplicated end knot must match master (:left == :right).
@testset "extrapolation past duplicated end knot unchanged" begin
    tend = 2.0
    cb = DiscreteCallback(
        (u, t, integrator) -> t == tend,
        integrator -> (integrator.u[1] += 1.0);
        save_positions = (true, true),
    )
    prob = ODEProblem((u, p, t) -> -u, [1.0], (0.0, tend))
    sol = solve(
        prob, Tsit5(); callback = cb, tstops = [tend], saveat = [0.0, tend],
        save_everystep = false,
    )
    @test sol.t[end - 1] == sol.t[end] == tend
    t_past = tend + 0.1
    @test sol(t_past; continuity = :left) == sol(t_past; continuity = :right)
    @test sol(t_past; continuity = :left) == sol.u[end]
end
