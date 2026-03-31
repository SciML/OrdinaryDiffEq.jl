using DiffEqBase: set_t!, set_u!, set_ut!
using OrdinaryDiffEq

println("First")
# set_X!(integrator, integrator.X) should not change the result.
@testset "Trivial $setter ($alg, inplace=$iip)" for alg in [RK4, Trapezoid],
        setter in [set_t!, set_u!, set_ut!],
        iip in [false, true]

    if iip
        f = (du, u, p, t) -> (du .= 2u)
        prob = ODEProblem{iip}(f, [1 / 2], (0.0, 2.0))
    else
        f = (u, p, t) -> (2u)
        prob = ODEProblem{iip}(f, 1 / 2, (0.0, 2.0))
    end
    t_half = 1.0

    integrator1 = init(prob, alg(); save_everystep = false)
    integrator2 = init(prob, alg(); save_everystep = false)

    step!(integrator1, t_half)  # "inexact" stepping w/o tstops
    if setter === set_t!
        set_t!(integrator1, integrator1.t)
    elseif setter === set_u!
        set_u!(integrator1, integrator1.u)
    else
        set_ut!(integrator1, integrator1.u, integrator1.t)
    end
    solve!(integrator1)

    solve!(integrator2)

    @test integrator1.t == integrator2.t
    if alg === Trapezoid
        rtol = integrator1.opts.reltol * 100
        atol = integrator1.opts.abstol
        @test integrator1.u ≈ integrator2.u rtol = rtol atol = atol
    else
        @test integrator1.u == integrator2.u
    end
end

println("Second")

@testset "Resolve with $setter ($alg, inplace=$iip)" for alg in [RK4, Trapezoid],
        setter in [set_t!, set_ut!],
        iip in [false, true]

    if iip
        f = (du, u, p, t) -> (du .= 2u * cos(2π * t))
        prob1 = ODEProblem{iip}(f, [1 / 2], (0.0, 1.0))
    else
        f = (u, p, t) -> (2u * cos(2π * t))
        prob1 = ODEProblem{iip}(f, 1 / 2, (0.0, 1.0))
    end
    prob2 = remake(prob1; tspan = (0.0, 2.0))
    integrator1 = init(prob1, alg(); save_everystep = false)
    integrator2 = init(prob2, alg(); save_everystep = false)

    solve!(integrator1)
    if setter === set_t!
        set_t!(integrator1, 0)
    else
        set_ut!(integrator1, integrator1.u, 0)
    end
    solve!(integrator1)

    solve!(integrator2)

    rtol = integrator1.opts.reltol
    atol = integrator1.opts.abstol
    if alg === Trapezoid
        rtol *= 100
    end
    @test integrator1.u ≈ integrator2.u rtol = rtol atol = atol
end

@testset "set_proposed_dt!" begin
    # 1) Setup problem
    f = (du, u, p, t) -> (du[1] = p)
    u0 = [1.0]
    tspan = (0.0, 5.0)
    p = -1
    prob = ODEProblem(f, u0, tspan, p)

    # 2) Setup callback
    condition(u, t, integrator) = u[1] <= 0.0
    # function that sets dt that works with adaptive solvers:
    function affect!(integrator)
        integrator.p = 1
        set_proposed_dt!(integrator, 0.5)
    end
    cb = DiscreteCallback(condition, affect!; save_positions = (true, true))

    sol = solve(prob, RK4(), dt = 0.1, adaptive = false, callback = cb)
    @test diff(sol.t)[1] == 0.1
    @test diff(sol.t)[13] == 0.5

    sol = solve(prob, RK4(), dt = 0.1, adaptive = true, callback = cb, dtmax = 1.0)
    @test diff(sol.t)[1] == 0.1
    @test 0.5 in diff(sol.t)
end
