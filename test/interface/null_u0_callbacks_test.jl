# Tests for solving ODEProblems with null u0 and callbacks
# https://github.com/SciML/ModelingToolkit.jl/issues/4078

using OrdinaryDiffEq, DiffEqCallbacks, Test

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
end
