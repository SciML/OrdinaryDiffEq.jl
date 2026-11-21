using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqTsit5
using Test

const prob = prob_dde_constant_1delay_scalar

@testset "Basic iterator" begin
    # compute the solution of the DDE
    sol = solve(
        prob, MethodOfSteps(BS3());
        dt = 1 // 2^(4), tstops = [1.5], saveat = 0.01, save_everystep = true
    )

    # initialize integrator
    integrator = init(
        prob, MethodOfSteps(BS3());
        dt = 1 // 2^(4), tstops = [1.5], saveat = 0.01, save_everystep = true
    )

    # perform one step
    step!(integrator)
    @test integrator.iter == 1

    # move to next grid point
    integrator.opts.advance_to_tstop = true
    step!(integrator)
    @test integrator.t == 1.5

    # solve the DDE
    solve!(integrator)
    @test integrator.t == 10
    @test integrator.sol(9) ≈ sol(9)

    # move to next grid point
    push!(integrator.opts.tstops, 15)
    step!(integrator)
    @test integrator.t == 15

    # move just one step
    integrator.opts.advance_to_tstop = false
    step!(integrator)
    @test integrator.t > 15
end

@testset "Advanced iterators" begin
    # move to grid point in one step
    integrator1 = init(
        prob, MethodOfSteps(Tsit5());
        dt = 1 // 2^(4), tstops = [0.5], advance_to_tstop = true
    )
    for integ in integrator1
        @test 0.5 ≤ integ.t ≤ 10
    end

    # move to grid point in one step and stop there
    integrator2 = init(
        prob, MethodOfSteps(Tsit5());
        dt = 1 // 2^(4), tstops = [0.5], advance_to_tstop = true,
        stop_at_next_tstop = true
    )
    for integ in integrator2
        @test integ.t == 0.5
    end
    integrator2([10; 20])

    # step and show intervals
    integrator3 = init(prob, MethodOfSteps(Tsit5()); dt = 1 // 2^(4), tstops = [0.5])
    for integ in integrator3
        @show integ.tprev, integ.t
    end
    integrator3([10; 20])

    # iterator for chosen time points
    integrator4 = init(prob, MethodOfSteps(Tsit5()); dt = 1 // 2^(4))
    ts = 1:10
    us = Float64[]
    for t in ts
        while integrator4.t < t && !isempty(integrator4.opts.tstops)
            step!(integrator4)
        end
        while integrator4.t < t
            step!(integrator4)
        end
        push!(us, integrator4(t))
    end
    @test us ≈ integrator4.sol(ts).u
end

@testset "iter" begin
    integrator = init(prob, MethodOfSteps(RK4()); dt = 1 // 2^(9), adaptive = false)
    for i in Base.Iterators.take(integrator, 12)
    end
    @test integrator.iter == 12
    for i in Base.Iterators.take(integrator, 12)
    end
    @test integrator.iter == 24
end
