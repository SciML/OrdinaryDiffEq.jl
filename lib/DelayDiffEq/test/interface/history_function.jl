using DelayDiffEq
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqCore
using Test

# check constant extrapolation with problem with vanishing delays at t = 0
@testset "vanishing delays" begin
    prob = DDEProblem((u, h, p, t) -> -h(p, t / 2), 1.0, (p, t) -> 1.0, (0.0, 10.0))
    solve(prob, MethodOfSteps(RK4()))
end

@testset "general" begin
    # naive history functions
    h_notinplace(p, t; idxs = nothing) = idxs === nothing ? [t, -t] : [t, -t][idxs]

    function h_inplace(val, p, t; idxs = nothing)
        if idxs === nothing
            val[1] = t
            val[2] = -t
        else
            val .= [t; -t][idxs]
        end
    end

    # ODE integrator
    prob = ODEProblem((du, u, p, t) -> @.(du = p * u), ones(2), (0.0, 1.0), 1.01)
    integrator = init(prob, Tsit5())

    # combined history function
    history_notinplace = DelayDiffEq.HistoryFunction(
        h_notinplace,
        integrator
    )
    history_inplace = DelayDiffEq.HistoryFunction(
        h_inplace,
        integrator
    )

    # test evaluation of history function
    @testset "evaluation" for idxs in (nothing, [2])
        # expected value
        trueval = h_notinplace(nothing, -1; idxs = idxs)

        # out-of-place
        @test history_notinplace(nothing, -1, Val{0}; idxs = idxs) == trueval

        # in-place
        val = zero(trueval)
        history_inplace(val, nothing, -1; idxs = idxs)
        @test val == trueval

        val = zero(trueval)
        history_inplace(val, nothing, -1, Val{0}; idxs = idxs)
        @test val == trueval
    end

    # test constant extrapolation
    @testset "constant extrapolation" for deriv in (Val{0}, Val{1}), idxs in (nothing, [2])
        # expected value
        trueval = deriv == Val{0} ?
            (idxs === nothing ? integrator.u : integrator.u[[2]]) :
            (idxs === nothing ? zeros(2) : [0.0])

        # out-of-place
        history_notinplace.isout = false
        @test history_notinplace(nothing, 1, deriv; idxs = idxs) == trueval
        @test history_notinplace.isout

        # in-place
        history_inplace.isout = false
        @test history_inplace(nothing, nothing, 1, deriv; idxs = idxs) == trueval
        @test history_inplace.isout

        history_inplace.isout = false
        val = 1 .- trueval # ensures that val ≠ trueval
        history_inplace(val, nothing, 1, deriv; idxs = idxs)
        @test val == trueval
        @test history_inplace.isout
    end

    # add step to integrator
    @testset "update integrator" begin
        OrdinaryDiffEqCore.loopheader!(integrator)
        OrdinaryDiffEqCore.perform_step!(integrator, integrator.cache)
        integrator.t = integrator.dt
        @test 0.01 < integrator.t < 1
        @test integrator.sol.t[end] == 0
    end

    # test integrator interpolation
    @testset "integrator interpolation" for deriv in (Val{0}, Val{1}),
            idxs in (nothing, [2])
        # expected value
        trueval = OrdinaryDiffEqCore.current_interpolant(0.01, integrator, idxs, deriv)

        # out-of-place
        history_notinplace.isout = false
        @test history_notinplace(nothing, 0.01, deriv; idxs = idxs) == trueval
        @test history_notinplace.isout

        # in-place
        history_inplace.isout = false
        val = zero(trueval)
        history_inplace(val, nothing, 0.01, deriv; idxs = idxs)
        @test val == trueval
        @test history_inplace.isout
    end

    # add step to solution
    @testset "update solution" begin
        integrator.t = 0
        OrdinaryDiffEqCore.loopfooter!(integrator)
        @test integrator.t == integrator.sol.t[end]
    end

    # test solution interpolation
    @testset "solution interpolation" for deriv in (Val{0}, Val{1}), idxs in (nothing, [2])
        # expected value
        trueval = integrator.sol.interp(0.01, idxs, deriv, integrator.p)

        # out-of-place
        history_notinplace.isout = false
        @test history_notinplace(nothing, 0.01, deriv; idxs = idxs) == trueval
        @test !history_notinplace.isout

        # in-place
        history_inplace.isout = false
        val = zero(trueval)
        history_inplace(val, nothing, 0.01, deriv; idxs = idxs)
        @test val == trueval
        @test !history_inplace.isout
    end

    # test integrator extrapolation
    @testset "integrator extrapolation" for deriv in (Val{0}, Val{1}), idxs in (0, [2])

        idxs == 0 && (idxs = nothing)
        # expected value
        trueval = OrdinaryDiffEqCore.current_interpolant(1, integrator, idxs, deriv)

        # out-of-place
        history_notinplace.isout = false
        @test history_notinplace(nothing, 1, deriv; idxs = idxs) == trueval
        @test history_notinplace.isout

        # in-place
        history_inplace.isout = false
        val = zero(trueval)
        history_inplace(val, nothing, 1, deriv; idxs = idxs)
        @test val == trueval
        @test history_inplace.isout
    end
end
