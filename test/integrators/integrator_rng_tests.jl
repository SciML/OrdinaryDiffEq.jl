using OrdinaryDiffEq, Test, Random, StableRNGs

# Simple ODE for testing: du/dt = 2u
f_oop(u, p, t) = 2u
f_iip(du, u, p, t) = (du .= 2 .* u)

@testset "Integrator RNG Interface" begin
    @testset "Default RNG (rng not provided)" begin
        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        integrator = init(prob, Tsit5())

        @test SciMLBase.has_rng(integrator)
        @test SciMLBase.get_rng(integrator) === Random.default_rng()
    end

    @testset "Custom RNG via init (out-of-place)" begin
        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        rng = Random.Xoshiro(42)
        integrator = init(prob, Tsit5(); rng)

        @test SciMLBase.has_rng(integrator)
        @test SciMLBase.get_rng(integrator) === rng
    end

    @testset "Custom RNG via init (in-place)" begin
        prob = ODEProblem(f_iip, [0.5], (0.0, 1.0))
        rng = Random.Xoshiro(123)
        integrator = init(prob, Tsit5(); rng)

        @test SciMLBase.has_rng(integrator)
        @test SciMLBase.get_rng(integrator) === rng
    end

    @testset "Custom RNG via solve propagates to integrator" begin
        rng_from_callback = Ref{Any}(nothing)
        cb = DiscreteCallback(
            (u, t, integrator) -> true,
            integrator -> begin
                rng_from_callback[] = SciMLBase.get_rng(integrator)
                return nothing
            end
        )

        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        rng = Random.Xoshiro(99)
        sol = solve(prob, Tsit5(); rng, callback = cb)

        @test sol.retcode == ReturnCode.Success
        @test rng_from_callback[] === rng
    end

    @testset "set_rng! replaces the RNG" begin
        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        rng1 = Random.Xoshiro(1)
        integrator = init(prob, Tsit5(); rng = rng1)

        @test SciMLBase.get_rng(integrator) === rng1

        rng2 = Random.Xoshiro(2)
        SciMLBase.set_rng!(integrator, rng2)
        @test SciMLBase.get_rng(integrator) === rng2
    end

    @testset "reinit! with rng kwarg" begin
        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        rng1 = Random.Xoshiro(10)
        integrator = init(prob, Tsit5(); rng = rng1)
        @test SciMLBase.get_rng(integrator) === rng1

        # reinit! without rng should keep existing RNG
        reinit!(integrator)
        @test SciMLBase.get_rng(integrator) === rng1

        # reinit! with rng should replace it
        rng2 = Random.Xoshiro(20)
        reinit!(integrator; rng = rng2)
        @test SciMLBase.get_rng(integrator) === rng2
    end

    @testset "reinit! with rng sets RNG before callback initialization" begin
        rng_seen_in_init = Ref{Any}(nothing)
        cb = DiscreteCallback(
            (u, t, integrator) -> true,
            integrator -> nothing;
            initialize = (cb, u, t, integrator) -> begin
                rng_seen_in_init[] = SciMLBase.get_rng(integrator)
                return nothing
            end
        )

        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        rng1 = Random.Xoshiro(10)
        integrator = init(prob, Tsit5(); rng = rng1, callback = cb)

        rng2 = Random.Xoshiro(20)
        reinit!(integrator; rng = rng2, reinit_callbacks = true)

        # The callback's initialize hook should see rng2, not rng1
        @test rng_seen_in_init[] === rng2
    end

    @testset "RNG preserved across solve! cycle" begin
        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        rng = Random.Xoshiro(42)
        integrator = init(prob, Tsit5(); rng)
        solve!(integrator)

        @test SciMLBase.get_rng(integrator) === rng
        @test integrator.sol.retcode == ReturnCode.Success
    end

    @testset "RNG type parameter is concrete" begin
        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        rng = Random.Xoshiro(42)
        integrator = init(prob, Tsit5(); rng)

        # The RNG type should be a concrete type parameter, not Any
        @test typeof(integrator).parameters[end] === typeof(rng)
    end

    @testset "Callback can access RNG via get_rng" begin
        rng_from_callback = Ref{Any}(nothing)
        cb = DiscreteCallback(
            (u, t, integrator) -> true,
            integrator -> begin
                rng_from_callback[] = SciMLBase.get_rng(integrator)
                return nothing
            end
        )

        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        rng = Random.Xoshiro(42)
        sol = solve(prob, Tsit5(); rng, callback = cb)

        @test rng_from_callback[] === rng
    end

    @testset "set_rng! with incompatible type throws ArgumentError" begin
        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        rng = Random.Xoshiro(42)
        integrator = init(prob, Tsit5(); rng)

        @test_throws ArgumentError SciMLBase.set_rng!(integrator, Random.MersenneTwister(1))
    end

    @testset "Different solvers support rng kwarg" begin
        prob_oop = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        prob_iip = ODEProblem(f_iip, [0.5], (0.0, 1.0))
        rng = Random.Xoshiro(42)

        for (alg, prob) in [
                (Tsit5(), prob_oop),
                (Vern7(), prob_oop),
                (RK4(), prob_iip),
                (Rosenbrock23(), prob_iip),
            ]
            integrator = init(prob, alg; rng)
            @test SciMLBase.has_rng(integrator)
            @test SciMLBase.get_rng(integrator) === rng
            solve!(integrator)
            @test integrator.sol.retcode == ReturnCode.Success
        end
    end

    @testset "Callback rand draws are reproducible with StableRNG" begin
        # A callback that draws from the integrator RNG at every step
        function make_collecting_callback()
            draws = Float64[]
            cb = DiscreteCallback(
                (u, t, integrator) -> true,
                integrator -> begin
                    push!(draws, rand(SciMLBase.get_rng(integrator)))
                    return nothing
                end
            )
            return cb, draws
        end

        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))

        cb1, draws1 = make_collecting_callback()
        solve(prob, Tsit5(); rng = StableRNG(42), callback = cb1)

        cb2, draws2 = make_collecting_callback()
        solve(prob, Tsit5(); rng = StableRNG(42), callback = cb2)

        @test !isempty(draws1)
        @test draws1 == draws2
    end

    @testset "Callback rand draws differ with different StableRNG seeds" begin
        function make_collecting_callback()
            draws = Float64[]
            cb = DiscreteCallback(
                (u, t, integrator) -> true,
                integrator -> begin
                    push!(draws, rand(SciMLBase.get_rng(integrator)))
                    return nothing
                end
            )
            return cb, draws
        end

        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))

        cb1, draws1 = make_collecting_callback()
        solve(prob, Tsit5(); rng = StableRNG(42), callback = cb1)

        cb2, draws2 = make_collecting_callback()
        solve(prob, Tsit5(); rng = StableRNG(99), callback = cb2)

        @test !isempty(draws1)
        @test draws1 != draws2
    end

    @testset "reinit! with new StableRNG resets rand sequence" begin
        draws = Float64[]
        cb = DiscreteCallback(
            (u, t, integrator) -> true,
            integrator -> begin
                push!(draws, rand(SciMLBase.get_rng(integrator)))
                return nothing
            end
        )

        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        integrator = init(prob, Tsit5(); rng = StableRNG(42), callback = cb)
        solve!(integrator)
        draws_run1 = copy(draws)

        empty!(draws)
        reinit!(integrator; rng = StableRNG(42))
        solve!(integrator)
        draws_run2 = draws

        @test !isempty(draws_run1)
        @test draws_run1 == draws_run2
    end

    @testset "StableRNG type parameter is concrete" begin
        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        rng = StableRNG(42)
        integrator = init(prob, Tsit5(); rng)

        @test typeof(integrator).parameters[end] === StableRNG
        @test SciMLBase.get_rng(integrator) === rng
    end

    @testset "set_rng! works with same StableRNG type" begin
        prob = ODEProblem(f_oop, 0.5, (0.0, 1.0))
        rng1 = StableRNG(1)
        integrator = init(prob, Tsit5(); rng = rng1)

        rng2 = StableRNG(2)
        SciMLBase.set_rng!(integrator, rng2)
        @test SciMLBase.get_rng(integrator) === rng2

        # Cross-type should fail
        @test_throws ArgumentError SciMLBase.set_rng!(integrator, Random.Xoshiro(1))
    end
end
