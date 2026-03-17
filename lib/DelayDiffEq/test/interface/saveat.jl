using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqTsit5
using Test

const prob = prob_dde_constant_1delay_long_ip
const alg = MethodOfSteps(Tsit5())

# reference integrator and solution
const dde_int = init(prob, alg)
const sol = solve!(dde_int)

@testset "reference" begin
    # solution equals solution of DDE integrator
    @test sol.t == dde_int.sol.t
    @test sol.u == dde_int.sol.u

    # solution equals solution of ODE integrator
    @test sol.t == dde_int.integrator.sol.t
    @test sol.u == dde_int.integrator.sol.u
end

# do not save every step
@testset "not every step (save_start=$save_start)" for save_start in (false, true)
    # for time(s) as scalar (implicitly adds end point as well!) and vectors
    for saveat in (25.0, [25.0, 50.0, 75.0])
        dde_int2 = init(prob, alg; saveat = saveat, save_start = save_start)

        # end point is saved if saveat is a scalar
        @test dde_int2.opts.save_end == (saveat isa Number)

        sol2 = solve!(dde_int2)

        # solution is equal to solution of DDE integrator
        @test sol2.t == dde_int2.sol.t
        @test sol2.u == dde_int2.sol.u

        # time point of solution
        if saveat isa Number
            @test sol2.t ==
                (save_start ? [0.0, 25.0, 50.0, 75.0, 100.0] : [25.0, 50.0, 75.0, 100.0])
        else
            @test sol2.t == (save_start ? [0.0, 25.0, 50.0, 75.0] : [25.0, 50.0, 75.0])
        end

        # history is equal to solution above
        @test sol.t == dde_int2.integrator.sol.t
        @test sol.u == dde_int2.integrator.sol.u
    end
end

# do not save every step
@testset "not every step (save_end=$save_end)" for save_end in (false, true)
    # for time(s) as scalar (implicitly adds end point as well!) and vectors
    for saveat in (25.0, [25.0, 50.0, 75.0])
        dde_int2 = init(prob, alg; saveat = saveat, save_end = save_end)

        # start point is saved if saveat is a scalar
        @test dde_int2.opts.save_start == (saveat isa Number)

        sol2 = solve!(dde_int2)

        # solution is equal to solution of DDE integrator
        @test sol2.t == dde_int2.sol.t
        @test sol2.u == dde_int2.sol.u

        # time point of solution
        if saveat isa Number
            @test sol2.t ==
                (save_end ? [0.0, 25.0, 50.0, 75.0, 100.0] : [0.0, 25.0, 50.0, 75.0])
        else
            @test sol2.t == (save_end ? [25.0, 50.0, 75.0, 100.0] : [25.0, 50.0, 75.0])
        end

        # history is equal to solution above
        @test sol.t == dde_int2.integrator.sol.t
        @test sol.u == dde_int2.integrator.sol.u
    end
end

# save every step
@testset "every step (save_start=$save_start)" for save_start in (false, true)
    for saveat in (25.0, [25.0, 50.0, 75.0])
        dde_int2 = init(
            prob, alg; saveat = saveat, save_everystep = true,
            save_start = save_start
        )

        # end point is saved implicitly
        @test dde_int2.opts.save_end

        sol2 = solve!(dde_int2)

        # solution is equal to solution of DDE integrator
        @test sol2.t == dde_int2.sol.t
        @test sol2.u == dde_int2.sol.u

        # time points of solution
        @test symdiff(sol.t, sol2.t) ==
            (save_start ? [25.0, 50.0, 75.0] : [0.0, 25.0, 50.0, 75.0])

        # history is equal to solution above
        @test sol.t == dde_int2.integrator.sol.t
        @test sol.u == dde_int2.integrator.sol.u
    end
end

# save every step
@testset "every step (save_end=$save_end)" for save_end in (false, true)
    for saveat in (25.0, [25.0, 50.0, 75.0])
        dde_int2 = init(
            prob, alg; saveat = saveat, save_everystep = true,
            save_end = save_end
        )

        # start point is saved implicitly
        @test dde_int2.opts.save_start

        sol2 = solve!(dde_int2)

        # solution is equal to solution of DDE integrator
        @test sol2.t == dde_int2.sol.t
        @test sol2.u == dde_int2.sol.u

        # time points of solution
        @test symdiff(sol.t, sol2.t) ==
            (save_end ? [25.0, 50.0, 75.0] : [100.0, 25.0, 50.0, 75.0])

        # history is equal to solution above
        @test sol.t == dde_int2.integrator.sol.t
        @test sol.u == dde_int2.integrator.sol.u
    end
end

@testset "not matching end time point" begin
    sol = solve(prob, alg; saveat = 40)

    @test sol.t == [0.0, 40, 80, 100]
end

@testset "changing end time point saveat" begin
    _saveat = [0.0, 0.25, 0.5, 1.0]
    integ = init(
        DDEProblem((u, h, p, t) -> u, 0.0, (p, t) -> 0.0, (0.0, 1.0)),
        MethodOfSteps(Tsit5()), saveat = _saveat
    )
    add_tstop!(integ, 2.0)
    solve!(integ)
    @test integ.sol.t == _saveat

    integ = init(
        DDEProblem((u, h, p, t) -> u, 0.0, (p, t) -> 0.0, (0.0, 1.0)),
        MethodOfSteps(Tsit5()), saveat = _saveat, save_end = true
    )
    add_tstop!(integ, 2.0)
    solve!(integ)
    @test integ.sol.t == [0.0, 0.25, 0.5, 1.0, 2.0]

    integ = init(
        DDEProblem((u, h, p, t) -> u, 0.0, (p, t) -> 0.0, (0.0, 1.0)),
        MethodOfSteps(Tsit5()), saveat = _saveat, save_end = false
    )
    add_tstop!(integ, 2.0)
    solve!(integ)
    @test integ.sol.t == [0.0, 0.25, 0.5]
end
