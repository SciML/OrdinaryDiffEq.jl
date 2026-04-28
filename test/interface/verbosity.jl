using OrdinaryDiffEqCore
using OrdinaryDiffEqCore: DEVerbosity
using OrdinaryDiffEq
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using OrdinaryDiffEqExtrapolation, OrdinaryDiffEqFIRK, OrdinaryDiffEqRosenbrock,
    OrdinaryDiffEqSDIRK
using ODEProblemLibrary: prob_ode_vanderpol_stiff
using Test
import OrdinaryDiffEqCore.SciMLLogging as SciMLLogging
using LinearSolve: LinearVerbosity
using NonlinearSolve: NonlinearVerbosity

@testset "DEVerbosity Tests" begin
    @testset "Default constructor" begin
        v1 = DEVerbosity()
        @test v1 isa DEVerbosity

        # Test solver verbosity
        @test v1.linear_verbosity == SciMLLogging.Minimal()
        @test v1.nonlinear_verbosity == SciMLLogging.Minimal()

        # Test error control group (default: WarnLevel except newton_convergence, step_rejected, step_accepted, convergence_limit)
        @test v1.dt_NaN == SciMLLogging.WarnLevel
        @test v1.init_NaN == SciMLLogging.WarnLevel
        @test v1.dense_output_saveat == SciMLLogging.WarnLevel
        @test v1.max_iters == SciMLLogging.WarnLevel
        @test v1.dt_min_unstable == SciMLLogging.WarnLevel
        @test v1.instability == SciMLLogging.WarnLevel
        @test v1.newton_convergence == SciMLLogging.Silent
        @test v1.step_rejected == SciMLLogging.Silent
        @test v1.step_accepted == SciMLLogging.Silent
        @test v1.convergence_limit == SciMLLogging.Silent

        # Test performance group (default: Silent except mismatched_input_output_type)
        @test v1.alg_switch == SciMLLogging.Silent
        @test v1.stiff_detection == SciMLLogging.Silent
        @test v1.mismatched_input_output_type == SciMLLogging.WarnLevel
        @test v1.jacobian_update == SciMLLogging.Silent
        @test v1.w_factorization == SciMLLogging.Silent
        @test v1.newton_iterations == SciMLLogging.Silent

        # Test numerical group (default: WarnLevel except shampine_dt, stability_check)
        @test v1.rosenbrock_no_differential_states == SciMLLogging.WarnLevel
        @test v1.shampine_dt == SciMLLogging.Silent
        @test v1.unlimited_dt == SciMLLogging.WarnLevel
        @test v1.dt_epsilon == SciMLLogging.WarnLevel
        @test v1.stability_check == SciMLLogging.Silent
        @test v1.near_singular == SciMLLogging.Silent
    end

    @testset "DEVerbosity preset constructors" begin
        v_none = DEVerbosity(SciMLLogging.None())
        v_minimal = DEVerbosity(SciMLLogging.Minimal())
        v_standard = DEVerbosity(SciMLLogging.Standard())
        v_detailed = DEVerbosity(SciMLLogging.Detailed())
        v_all = DEVerbosity(SciMLLogging.All())

        # Test None - everything Silent
        @test v_none.linear_verbosity == SciMLLogging.None()
        @test v_none.nonlinear_verbosity == SciMLLogging.None()
        @test v_none.dt_NaN == SciMLLogging.Silent
        @test v_none.init_NaN == SciMLLogging.Silent
        @test v_none.alg_switch == SciMLLogging.Silent
        @test v_none.rosenbrock_no_differential_states == SciMLLogging.Silent

        # Test Minimal - only critical errors
        @test v_minimal.linear_verbosity == SciMLLogging.Minimal()
        @test v_minimal.nonlinear_verbosity == SciMLLogging.Minimal()
        @test v_minimal.dt_NaN == SciMLLogging.WarnLevel
        @test v_minimal.init_NaN == SciMLLogging.WarnLevel
        @test v_minimal.max_iters == SciMLLogging.WarnLevel
        @test v_minimal.newton_convergence == SciMLLogging.WarnLevel
        @test v_minimal.near_singular == SciMLLogging.WarnLevel
        @test v_minimal.alg_switch == SciMLLogging.Silent
        @test v_minimal.dense_output_saveat == SciMLLogging.Silent
        @test v_minimal.mismatched_input_output_type == SciMLLogging.Silent

        # Test Standard - same as default
        @test v_standard.dt_NaN == SciMLLogging.WarnLevel
        @test v_standard.init_NaN == SciMLLogging.WarnLevel
        @test v_standard.alg_switch == SciMLLogging.Silent
        @test v_standard.dense_output_saveat == SciMLLogging.WarnLevel
        @test v_standard.mismatched_input_output_type == SciMLLogging.WarnLevel

        # Test Detailed - includes debugging info
        @test v_detailed.linear_verbosity == SciMLLogging.Detailed()
        @test v_detailed.nonlinear_verbosity == SciMLLogging.Detailed()
        @test v_detailed.alg_switch == SciMLLogging.InfoLevel
        @test v_detailed.dense_output_saveat == SciMLLogging.InfoLevel
        @test v_detailed.shampine_dt == SciMLLogging.InfoLevel
        @test v_detailed.jacobian_update == SciMLLogging.InfoLevel
        @test v_detailed.w_factorization == SciMLLogging.InfoLevel
        @test v_detailed.convergence_limit == SciMLLogging.InfoLevel

        # Test All - maximum verbosity
        @test v_all.linear_verbosity == SciMLLogging.All()
        @test v_all.nonlinear_verbosity == SciMLLogging.All()
        @test v_all.alg_switch == SciMLLogging.InfoLevel
        @test v_all.shampine_dt == SciMLLogging.InfoLevel
        @test v_all.dense_output_saveat == SciMLLogging.InfoLevel
        @test v_all.step_rejected == SciMLLogging.InfoLevel
        @test v_all.step_accepted == SciMLLogging.InfoLevel
        @test v_all.stiff_detection == SciMLLogging.InfoLevel
    end

    @testset "Group-level keyword constructors" begin
        v_error = DEVerbosity(error_control = SciMLLogging.ErrorLevel)
        # Test all error_control fields
        @test v_error.dt_NaN == SciMLLogging.ErrorLevel
        @test v_error.init_NaN == SciMLLogging.ErrorLevel
        @test v_error.dense_output_saveat == SciMLLogging.ErrorLevel
        @test v_error.max_iters == SciMLLogging.ErrorLevel
        @test v_error.dt_min_unstable == SciMLLogging.ErrorLevel
        @test v_error.instability == SciMLLogging.ErrorLevel
        @test v_error.newton_convergence == SciMLLogging.ErrorLevel
        @test v_error.step_rejected == SciMLLogging.ErrorLevel
        @test v_error.step_accepted == SciMLLogging.ErrorLevel
        @test v_error.convergence_limit == SciMLLogging.ErrorLevel

        v_numerical = DEVerbosity(numerical = SciMLLogging.Silent)
        # Test all numerical fields
        @test v_numerical.rosenbrock_no_differential_states == SciMLLogging.Silent
        @test v_numerical.shampine_dt == SciMLLogging.Silent
        @test v_numerical.unlimited_dt == SciMLLogging.Silent
        @test v_numerical.dt_epsilon == SciMLLogging.Silent
        @test v_numerical.stability_check == SciMLLogging.Silent
        @test v_numerical.near_singular == SciMLLogging.Silent

        v_performance = DEVerbosity(performance = SciMLLogging.InfoLevel)
        # Test all performance fields
        @test v_performance.alg_switch == SciMLLogging.InfoLevel
        @test v_performance.stiff_detection == SciMLLogging.InfoLevel
        @test v_performance.mismatched_input_output_type == SciMLLogging.InfoLevel
        @test v_performance.jacobian_update == SciMLLogging.InfoLevel
        @test v_performance.w_factorization == SciMLLogging.InfoLevel
        @test v_performance.newton_iterations == SciMLLogging.InfoLevel
    end

    @testset "Mixed group and individual settings" begin
        v_mixed = DEVerbosity(
            numerical = SciMLLogging.Silent,
            shampine_dt = SciMLLogging.WarnLevel,
            performance = SciMLLogging.InfoLevel
        )
        # Individual override should take precedence
        @test v_mixed.shampine_dt == SciMLLogging.WarnLevel
        # Other numerical options should use group setting
        @test v_mixed.rosenbrock_no_differential_states == SciMLLogging.Silent
        @test v_mixed.unlimited_dt == SciMLLogging.Silent
        # Performance group setting should apply
        @test v_mixed.alg_switch == SciMLLogging.InfoLevel
        @test v_mixed.mismatched_input_output_type == SciMLLogging.InfoLevel
    end

    @testset "Individual keyword arguments" begin
        v_individual = DEVerbosity(
            dt_NaN = SciMLLogging.ErrorLevel,
            alg_switch = SciMLLogging.InfoLevel,
            shampine_dt = SciMLLogging.Silent
        )
        @test v_individual.dt_NaN == SciMLLogging.ErrorLevel
        @test v_individual.alg_switch == SciMLLogging.InfoLevel
        @test v_individual.shampine_dt == SciMLLogging.Silent
        # Unspecified options should use defaults
        @test v_individual.init_NaN == SciMLLogging.WarnLevel
        @test v_individual.unlimited_dt == SciMLLogging.WarnLevel
    end

    @testset "Linear and nonlinear verbosity passthrough" begin
        v_with_solvers = DEVerbosity(
            linear_verbosity = SciMLLogging.Detailed(),
            nonlinear_verbosity = SciMLLogging.Minimal()
        )
        @test v_with_solvers.linear_verbosity == SciMLLogging.Detailed()
        @test v_with_solvers.nonlinear_verbosity == SciMLLogging.Minimal()

        v_with_solvers2 = DEVerbosity(
            linear_verbosity = SciMLLogging.None(),
            nonlinear_verbosity = SciMLLogging.All()
        )
        @test v_with_solvers2.linear_verbosity == SciMLLogging.None()
        @test v_with_solvers2.nonlinear_verbosity == SciMLLogging.All()
    end

    @testset "Validation tests" begin
        # Test that invalid group arguments throw errors
        @test_throws ArgumentError DEVerbosity(error_control = "invalid")
        @test_throws ArgumentError DEVerbosity(performance = 123)
        @test_throws ArgumentError DEVerbosity(numerical = :wrong)

        # Test that invalid individual fields throw errors
        @test_throws ArgumentError DEVerbosity(dt_NaN = "invalid")
        @test_throws ArgumentError DEVerbosity(unknown_field = SciMLLogging.InfoLevel())

        # Test that Bool verbose is no longer supported (v7 breaking change)
        prob_simple = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
        @test_throws ArgumentError solve(prob_simple, Tsit5(); verbose = true)
        @test_throws ArgumentError solve(prob_simple, Tsit5(); verbose = false)
    end

    @testset "Multiple group settings" begin
        v = DEVerbosity(
            error_control = SciMLLogging.ErrorLevel,
            performance = SciMLLogging.InfoLevel,
            numerical = SciMLLogging.Silent
        )
        @test v.dt_NaN == SciMLLogging.ErrorLevel
        @test v.alg_switch == SciMLLogging.InfoLevel
        @test v.shampine_dt == SciMLLogging.Silent
    end

    @testset "Complex mixed settings" begin
        v = DEVerbosity(
            error_control = SciMLLogging.WarnLevel,
            performance = SciMLLogging.InfoLevel,
            numerical = SciMLLogging.Silent,
            linear_verbosity = SciMLLogging.Detailed(),
            nonlinear_verbosity = SciMLLogging.Minimal(),
            dt_NaN = SciMLLogging.ErrorLevel,  # Override specific error_control field
            shampine_dt = SciMLLogging.WarnLevel  # Override specific numerical field
        )
        # Check overrides took precedence
        @test v.dt_NaN == SciMLLogging.ErrorLevel
        @test v.shampine_dt == SciMLLogging.WarnLevel
        # Check other fields follow group settings
        @test v.init_NaN == SciMLLogging.WarnLevel
        @test v.alg_switch == SciMLLogging.InfoLevel
        @test v.unlimited_dt == SciMLLogging.Silent
        # Check solver verbosity
        @test v.linear_verbosity == SciMLLogging.Detailed()
        @test v.nonlinear_verbosity == SciMLLogging.Minimal()
    end

    @testset "Stiff Switching Message" begin
        verb = DEVerbosity(alg_switch = SciMLLogging.InfoLevel)
        solve(prob_ode_vanderpol_stiff, AutoTsit5(Rodas5()), verbose = verb)
    end

    @testset "Linear Verbosity Passthrough to Caches" begin
        # Define a simple stiff test problem
        function rober(du, u, p, t)
            y₁, y₂, y₃ = u
            k₁, k₂, k₃ = p
            du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
            du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
            du[3] = y₁ + y₂ + y₃ - 1
            nothing
        end
        u0 = [1.0, 0.0, 0.0]
        tspan = (0.0, 1.0e-1)
        p = [0.04, 3.0e7, 1.0e4]
        prob = ODEProblem(rober, u0, tspan, p)

        @testset "Rosenbrock Solvers" begin
            @testset "Rosenbrock23 with Detailed LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.Detailed())
                integrator = init(prob, Rosenbrock23(), verbose = verbose, dt = 1.0e-3)

                # Check that the cache has a linsolve field
                @test hasproperty(integrator.cache, :linsolve)

                # Check that the linear solver cache has verbose field
                @test hasproperty(integrator.cache.linsolve, :verbose)

                # Verify the verbosity was passed through correctly
                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Detailed())
            end

            @testset "Rosenbrock23 with Minimal LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.Minimal())
                integrator = init(prob, Rosenbrock23(), verbose = verbose, dt = 1.0e-3)

                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Minimal())
            end

            @testset "Rosenbrock23 with None LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.None())
                integrator = init(prob, Rosenbrock23(), verbose = verbose, dt = 1.0e-3)

                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.None())
            end

            @testset "Rosenbrock32 with Detailed LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.Detailed())
                integrator = init(prob, Rosenbrock32(), verbose = verbose, dt = 1.0e-3)

                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Detailed())
            end

            @testset "Rodas4 with All LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.All())
                integrator = init(prob, Rodas4(), verbose = verbose, dt = 1.0e-3)

                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.All())
            end
        end

        @testset "FIRK Solvers (Radau Methods)" begin
            @testset "RadauIIA3 with Detailed LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.Detailed())
                integrator = init(prob, RadauIIA3(), verbose = verbose, dt = 1.0e-3)

                # RadauIIA3 has a linsolve field
                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Detailed())
            end

            @testset "RadauIIA3 with Minimal LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.Minimal())
                integrator = init(prob, RadauIIA3(), verbose = verbose, dt = 1.0e-3)

                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Minimal())
            end

            @testset "RadauIIA5 with Detailed LinearVerbosity (two linear solvers)" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.Detailed())
                integrator = init(prob, RadauIIA5(), verbose = verbose, dt = 1.0e-3)

                @test integrator.cache.linsolve1.verbose == LinearVerbosity(SciMLLogging.Detailed())
                @test integrator.cache.linsolve2.verbose == LinearVerbosity(SciMLLogging.Detailed())
            end

            @testset "RadauIIA5 with None LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.None())
                integrator = init(prob, RadauIIA5(), verbose = verbose, dt = 1.0e-3)

                @test integrator.cache.linsolve1.verbose == LinearVerbosity(SciMLLogging.None())
                @test integrator.cache.linsolve2.verbose == LinearVerbosity(SciMLLogging.None())
            end

            @testset "RadauIIA9 with All LinearVerbosity (three linear solvers)" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.All())
                integrator = init(prob, RadauIIA9(), verbose = verbose, dt = 1.0e-3)

                # Check all three linear solvers have the correct verbosity
                @test integrator.cache.linsolve1.verbose == LinearVerbosity(SciMLLogging.All())
                @test integrator.cache.linsolve2.verbose ==
                    LinearVerbosity(SciMLLogging.All())
                @test integrator.cache.linsolve3.verbose ==
                    LinearVerbosity(SciMLLogging.All())
            end
        end

        @testset "Extrapolation Solvers" begin
            @testset "ImplicitEulerExtrapolation with Detailed LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.Detailed())
                integrator = init(prob, ImplicitEulerExtrapolation(), verbose = verbose, dt = 1.0e-3)

                # Check that all linear solvers in the array have the correct verbosity
                for ls in integrator.cache.linsolve
                    @test ls.verbose == LinearVerbosity(SciMLLogging.Detailed())
                end
            end

            @testset "ImplicitEulerExtrapolation with Minimal LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.Minimal())
                integrator = init(prob, ImplicitEulerExtrapolation(), verbose = verbose, dt = 1.0e-3)

                for ls in integrator.cache.linsolve
                    @test ls.verbose == LinearVerbosity(SciMLLogging.Minimal())
                end
            end

            @testset "ImplicitDeuflhardExtrapolation with None LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.None())
                integrator = init(prob, ImplicitDeuflhardExtrapolation(), verbose = verbose, dt = 1.0e-3)

                for ls in integrator.cache.linsolve
                    @test ls.verbose == LinearVerbosity(SciMLLogging.None())
                end
            end

            @testset "ImplicitHairerWannerExtrapolation with All LinearVerbosity" begin
                verbose = DEVerbosity(linear_verbosity = SciMLLogging.All())
                integrator = init(prob, ImplicitHairerWannerExtrapolation(), verbose = verbose, dt = 1.0e-3)

                for ls in integrator.cache.linsolve
                    @test ls.verbose == LinearVerbosity(SciMLLogging.All())
                end
            end
        end

        @testset "Preset Verbosity Levels" begin
            @testset "Rosenbrock23 with Standard() preset (default linear_verbosity = Minimal)" begin
                verbose = DEVerbosity(SciMLLogging.Standard())
                integrator = init(prob, Rosenbrock23(), verbose = verbose, dt = 1.0e-3)

                # Standard() uses Minimal() for linear_verbosity
                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Minimal())
            end

            @testset "RadauIIA3 with Detailed() preset" begin
                verbose = DEVerbosity(SciMLLogging.Detailed())
                integrator = init(prob, RadauIIA3(), verbose = verbose, dt = 1.0e-3)

                # Detailed() uses Detailed() for linear_verbosity
                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Detailed())
            end

        end
    end

    @testset "Nonlinear Verbosity Passthrough to Caches" begin
        # Define a simple stiff test problem
        function rober(du, u, p, t)
            y₁, y₂, y₃ = u
            k₁, k₂, k₃ = p
            du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
            du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
            du[3] = y₁ + y₂ + y₃ - 1
            nothing
        end
        u0 = [1.0, 0.0, 0.0]
        tspan = (0.0, 1.0e-1)
        p = [0.04, 3.0e7, 1.0e4]
        prob = ODEProblem(rober, u0, tspan, p)

        @testset "ImplicitEuler with Detailed NonlinearVerbosity" begin
            verbose = DEVerbosity(nonlinear_verbosity = SciMLLogging.Detailed())
            integrator = init(prob, ImplicitEuler(nlsolve = NonlinearSolveAlg()), verbose = verbose, dt = 1.0e-3)

            # Check that the cache has an nlsolver field
            @test hasproperty(integrator.cache, :nlsolver)

            # Check that the nlsolver has a cache field
            @test hasproperty(integrator.cache.nlsolver, :cache)

            # Verify the verbosity was passed through correctly
            @test integrator.cache.nlsolver.cache.cache.verbose == NonlinearVerbosity(SciMLLogging.Detailed())
        end

    end
end
