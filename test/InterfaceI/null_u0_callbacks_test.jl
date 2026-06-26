# Tests for solving ODEProblems with null u0 and callbacks
# https://github.com/SciML/ModelingToolkit.jl/issues/4078

using OrdinaryDiffEq, DiffEqCallbacks, Test
using OrdinaryDiffEqFunctionMap, OrdinaryDiffEqTsit5, OrdinaryDiffEqLowOrderRK,
    OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF

@testset "Null u0 with callbacks" begin
    @testset "Explicit solvers" begin
        for (name, alg, kwargs) in [
                ("FunctionMap", FunctionMap(), (dt = 0.1,)),
                ("Euler", Euler(), (dt = 0.1,)),
                ("RK4", RK4(), (dt = 0.1,)),
                ("Tsit5", Tsit5(), ()),
            ]
            @testset "$name" begin
                counter = Ref(0)
                cb = PresetTimeCallback([0.5], integrator -> (counter[] += 1))
                prob = ODEProblem(Returns(nothing), nothing, (0.0, 1.0))

                sol = solve(prob, alg; callback = cb, kwargs...)

                @test sol.retcode == ReturnCode.Success
                @test counter[] == 1  # Callback triggered once at t=0.5
                @test 0.5 in sol.t    # Solution includes callback time
            end
        end
    end

    @testset "Implicit solvers" begin
        for (name, alg) in [
                ("Rosenbrock23", Rosenbrock23()),
                ("Rodas5P", Rodas5P()),
                ("TRBDF2", TRBDF2()),
                ("ImplicitEuler", ImplicitEuler()),
                ("QNDF", QNDF()),
            ]
            @testset "$name" begin
                counter = Ref(0)
                cb = PresetTimeCallback([0.5], integrator -> (counter[] += 1))
                prob = ODEProblem(Returns(nothing), nothing, (0.0, 1.0))

                sol = solve(prob, alg; callback = cb)

                @test sol.retcode == ReturnCode.Success
                @test counter[] == 1  # Callback triggered once at t=0.5
                @test 0.5 in sol.t    # Solution includes callback time
            end
        end
    end

    @testset "Multiple callbacks" begin
        counter = Ref(0)
        cb = PresetTimeCallback([0.25, 0.5, 0.75], integrator -> (counter[] += 1))
        prob = ODEProblem(Returns(nothing), nothing, (0.0, 1.0))

        sol = solve(prob, Tsit5(); callback = cb)

        @test sol.retcode == ReturnCode.Success
        @test counter[] == 3
    end

    @testset "DiscreteCallback with condition" begin
        counter = Ref(0)
        cb = DiscreteCallback(
            (u, t, integrator) -> t >= 0.5 && counter[] == 0,
            integrator -> (counter[] += 1)
        )
        prob = ODEProblem(Returns(nothing), nothing, (0.0, 1.0))

        sol = solve(prob, Tsit5(); callback = cb)

        @test sol.retcode == ReturnCode.Success
        @test counter[] == 1
    end

    @testset "Solution structure with null u0" begin
        prob = ODEProblem(Returns(nothing), nothing, (0.0, 1.0))
        sol = solve(prob, Tsit5())

        @test sol.retcode == ReturnCode.Success
        @test all(u -> u == Float64[], sol.u)  # All states are empty arrays
        @test sol.t[1] == 0.0
        @test sol.t[end] == 1.0
    end

    # https://github.com/SciML/OrdinaryDiffEq.jl/issues/3778
    @testset "Interpolating null u0" begin
        f!(du, u, p, t) = nothing
        prob = ODEProblem(f!, Float64[], (0.0, 1.0))
        sol = solve(prob, Tsit5())

        @test sol.retcode == ReturnCode.Success
        @test sol(0.5; idxs = Int[]) == Float64[]
        @test sol(0.5) == Float64[]
        @test sol(0.5) isa Vector{Float64}
        @test sol(0.5, Val{1}) == Float64[]  # derivative path
        # vector of times
        @test all(==(Float64[]), sol([0.25, 0.5, 0.75]).u)

        # in-step interpolation through the integrator hits the same path
        integrator = init(prob, Tsit5())
        step!(integrator)
        @test integrator(integrator.t - 1.0e-3) == Float64[]
    end
end
