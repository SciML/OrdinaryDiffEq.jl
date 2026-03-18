using StochasticDiffEq, Test, Random
using SDEProblemLibrary: prob_sde_linear

# Callable tstops support (e.g. SymbolicTstops from ModelingToolkit)
# The callable has signature (p, tspan) -> collection_of_tstop_times

@testset "Callable tstops" begin

    @testset "Basic callable tstops with EM (fixed step)" begin
        prob = prob_sde_linear
        my_tstops = (p, tspan) -> [0.33, 0.66]
        sol = solve(prob, EM(); dt = 1 // 2^4, tstops = my_tstops)
        @test 0.33 ∈ sol.t
        @test 0.66 ∈ sol.t
        @test sol.t[end] == prob.tspan[end]
    end

    @testset "Callable tstops with adaptive solver (SRIW1)" begin
        prob = prob_sde_linear
        my_tstops = (p, tspan) -> [0.33, 0.77]
        sol = solve(prob, SRIW1(); tstops = my_tstops)
        @test 0.33 ∈ sol.t
        @test 0.77 ∈ sol.t
        @test sol.t[end] == prob.tspan[end]
    end

    @testset "Callable tstops with ImplicitEM" begin
        f(u, p, t) = -p[1] * u
        g(u, p, t) = 0.1 * u
        prob = SDEProblem(f, g, 1.0, (0.0, 10.0), [0.5])
        my_tstops = (p, tspan) -> [3.0, 6.0]
        sol = solve(prob, ImplicitEM(); dt = 0.01, tstops = my_tstops)
        @test 3.0 ∈ sol.t
        @test 6.0 ∈ sol.t
        @test sol.t[end] == 10.0
    end

    @testset "Callable tstops depending on parameters" begin
        f(u, p, t) = -p[1] * u
        g(u, p, t) = p[2] * u
        prob = SDEProblem(f, g, 1.0, (0.0, 5.0), [0.5, 0.1])

        # tstops that depend on parameter values
        my_tstops = (p, tspan) -> [p[1] * 2.0, p[1] * 4.0]
        sol = solve(prob, EM(); dt = 0.01, tstops = my_tstops)
        @test 1.0 ∈ sol.t  # 0.5 * 2.0
        @test 2.0 ∈ sol.t  # 0.5 * 4.0
        @test sol.t[end] == 5.0
    end

    @testset "Callable tstops depending on tspan" begin
        prob = prob_sde_linear
        # tstops at midpoint of tspan
        my_tstops = (p, tspan) -> [(tspan[1] + tspan[2]) / 2]
        sol = solve(prob, EM(); dt = 1 // 2^4, tstops = my_tstops)
        @test 0.5 ∈ sol.t
    end

    @testset "Callable tstops returning empty collection" begin
        prob = prob_sde_linear
        my_tstops = (p, tspan) -> Float64[]
        sol = solve(prob, EM(); dt = 1 // 2^4, tstops = my_tstops)
        @test sol.t[end] == prob.tspan[end]
    end

    @testset "Callable tstops with only tstops (no dt, fixed step)" begin
        # Fixed-step EM with callable tstops still requires dt
        prob = prob_sde_linear
        my_tstops = (p, tspan) -> [0.33, 0.66, 1.0]
        sol = solve(prob, EM(); dt = 1 // 2^4, tstops = my_tstops)
        @test 0.33 ∈ sol.t
        @test 0.66 ∈ sol.t
        @test 1.0 ∈ sol.t
    end

    @testset "Callable tstops with adaptive solver (no dt)" begin
        prob = prob_sde_linear
        my_tstops = (p, tspan) -> [0.33, 0.77]
        sol = solve(prob, SRIW1(); tstops = my_tstops)
        @test 0.33 ∈ sol.t
        @test 0.77 ∈ sol.t
    end

    @testset "Callable tstops with init/solve! workflow" begin
        prob = prob_sde_linear
        my_tstops = (p, tspan) -> [0.33, 0.77]
        integrator = init(prob, EM(); dt = 1 // 2^4, tstops = my_tstops)
        solve!(integrator)
        @test 0.33 ∈ integrator.sol.t
        @test 0.77 ∈ integrator.sol.t
    end

    @testset "Callable tstops are stored in tstops_cache" begin
        prob = prob_sde_linear
        my_tstops = (p, tspan) -> [0.33]
        integrator = init(prob, EM(); dt = 1 // 2^4, tstops = my_tstops)
        # The callable should be stored in tstops_cache for reinit!
        @test integrator.opts.tstops_cache === my_tstops
    end

    @testset "Callable tstops with reinit!" begin
        f(u, p, t) = -p[1] * u
        g(u, p, t) = 0.1 * u
        prob = SDEProblem(f, g, 1.0, (0.0, 5.0), [0.5])

        # Parameter-dependent callable tstops
        my_tstops = (p, tspan) -> [p[1] * 4.0]
        integrator = init(prob, EM(); dt = 0.01, tstops = my_tstops)
        solve!(integrator)
        @test 2.0 ∈ integrator.sol.t  # 0.5 * 4.0

        # Reinit reuses the callable from tstops_cache
        reinit!(integrator)
        solve!(integrator)
        @test 2.0 ∈ integrator.sol.t  # 0.5 * 4.0 again
    end

    @testset "Callable tstops with reinit! picks up changed parameters" begin
        f(u, p, t) = -p[1] * u
        g(u, p, t) = 0.1 * u
        prob = SDEProblem(f, g, 1.0, (0.0, 5.0), [0.5])

        # tstop at p[1] * 4.0
        my_tstops = (p, tspan) -> [p[1] * 4.0]
        integrator = init(prob, EM(); dt = 0.01, tstops = my_tstops)
        solve!(integrator)
        @test 2.0 ∈ integrator.sol.t   # 0.5 * 4.0
        @test 3.0 ∉ integrator.sol.t

        # Change parameter before reinit!
        integrator.p[1] = 0.75
        reinit!(integrator)
        solve!(integrator)
        @test 3.0 ∈ integrator.sol.t   # 0.75 * 4.0
        @test 2.0 ∉ integrator.sol.t
    end

    @testset "allows_late_binding_tstops trait" begin
        @test SciMLBase.allows_late_binding_tstops(EM()) == true
        @test SciMLBase.allows_late_binding_tstops(SRIW1()) == true
        @test SciMLBase.allows_late_binding_tstops(ImplicitEM()) == true
        @test SciMLBase.allows_late_binding_tstops(LambaEM()) == true
        @test SciMLBase.allows_late_binding_tstops(EulerHeun()) == true
    end

    @testset "Callable tstops with in-place SDE" begin
        function f!(du, u, p, t)
            du[1] = -p[1] * u[1]
        end
        function g!(du, u, p, t)
            du[1] = 0.1 * u[1]
        end
        prob = SDEProblem(f!, g!, [1.0], (0.0, 5.0), [0.5])
        my_tstops = (p, tspan) -> [1.5, 3.0]
        sol = solve(prob, EM(); dt = 0.01, tstops = my_tstops)
        @test 1.5 ∈ sol.t
        @test 3.0 ∈ sol.t
    end

    @testset "Callable tstops with RODE solver" begin
        f(u, p, t, W) = 1.01u + 0.87u * W
        prob = RODEProblem(f, 1.0, (0.0, 5.0))
        my_tstops = (p, tspan) -> [2.0, 4.0]
        sol = solve(prob, RandomEM(); dt = 0.01, tstops = my_tstops)
        @test 2.0 ∈ sol.t
        @test 4.0 ∈ sol.t
    end

    @testset "Callable tstops with callbacks" begin
        # Zero drift and noise so u stays at its initial/reset value
        f(u, p, t) = 0.0
        g(u, p, t) = 0.0
        prob = SDEProblem(f, g, 1.0, (0.0, 5.0))

        # Discrete callback at t=2.0 that sets u to exactly 42.0
        condition(u, t, integrator) = t == 2.0
        affect!(integrator) = (integrator.u = 42.0)
        cb = DiscreteCallback(condition, affect!)

        my_tstops = (p, tspan) -> [2.0, 4.0]
        sol = solve(prob, EM(); dt = 0.01, tstops = my_tstops, callback = cb)
        @test 2.0 ∈ sol.t
        @test 4.0 ∈ sol.t
        # Verify callback fired: u should be 1.0 before t=2.0 and 42.0 after
        idx_before = findfirst(==(2.0), sol.t) - 1
        idx_after = findfirst(==(2.0), sol.t) + 1
        @test sol.u[idx_before] == 1.0
        @test sol.u[idx_after] == 42.0
        @test sol.u[end] == 42.0
    end

    @testset "Callable tstops with multiple solvers" begin
        f(u, p, t) = -0.5 * u
        g(u, p, t) = 0.1 * u
        prob = SDEProblem(f, g, 1.0, (0.0, 5.0))
        my_tstops = (p, tspan) -> [1.0, 3.0]
        for alg in (EM(), LambaEM(), EulerHeun())
            sol = solve(prob, alg; dt = 0.01, tstops = my_tstops)
            @test 1.0 ∈ sol.t
            @test 3.0 ∈ sol.t
        end
    end
end
