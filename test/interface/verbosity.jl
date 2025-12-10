using OrdinaryDiffEqCore
using OrdinaryDiffEqCore: ODEVerbosity
using OrdinaryDiffEq
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using ODEProblemLibrary: prob_ode_vanderpol_stiff
using Test
import OrdinaryDiffEqCore.SciMLLogging as SciMLLogging
using LinearSolve: LinearVerbosity
using NonlinearSolve: NonlinearVerbosity

@testset "ODEVerbosity Tests" begin
    @testset "Default constructor" begin
        v1 = ODEVerbosity()
        @test v1 isa ODEVerbosity

        # Test solver verbosity
        @test v1.linear_verbosity isa SciMLLogging.Minimal
        @test v1.nonlinear_verbosity isa SciMLLogging.Minimal

        # Test error control group (default: WarnLevel except newton_convergence, step_rejected, step_accepted, convergence_limit)
        @test v1.dt_NaN isa SciMLLogging.WarnLevel
        @test v1.init_NaN isa SciMLLogging.WarnLevel
        @test v1.dense_output_saveat isa SciMLLogging.WarnLevel
        @test v1.max_iters isa SciMLLogging.WarnLevel
        @test v1.dt_min_unstable isa SciMLLogging.WarnLevel
        @test v1.instability isa SciMLLogging.WarnLevel
        @test v1.newton_convergence isa SciMLLogging.Silent
        @test v1.step_rejected isa SciMLLogging.Silent
        @test v1.step_accepted isa SciMLLogging.Silent
        @test v1.convergence_limit isa SciMLLogging.Silent

        # Test performance group (default: Silent except mismatched_input_output_type)
        @test v1.alg_switch isa SciMLLogging.Silent
        @test v1.stiff_detection isa SciMLLogging.Silent
        @test v1.mismatched_input_output_type isa SciMLLogging.WarnLevel
        @test v1.jacobian_update isa SciMLLogging.Silent
        @test v1.w_factorization isa SciMLLogging.Silent
        @test v1.newton_iterations isa SciMLLogging.Silent

        # Test numerical group (default: WarnLevel except shampine_dt, dt_epsilon, stability_check)
        @test v1.rosenbrock_no_differential_states isa SciMLLogging.WarnLevel
        @test v1.shampine_dt isa SciMLLogging.Silent
        @test v1.unlimited_dt isa SciMLLogging.WarnLevel
        @test v1.dt_epsilon isa SciMLLogging.Silent
        @test v1.stability_check isa SciMLLogging.Silent
        @test v1.near_singular isa SciMLLogging.Silent
    end

    @testset "ODEVerbosity preset constructors" begin
        v_none = ODEVerbosity(SciMLLogging.None())
        v_minimal = ODEVerbosity(SciMLLogging.Minimal())
        v_standard = ODEVerbosity(SciMLLogging.Standard())
        v_detailed = ODEVerbosity(SciMLLogging.Detailed())
        v_all = ODEVerbosity(SciMLLogging.All())

        # Test None - everything Silent
        @test v_none.linear_verbosity isa SciMLLogging.None
        @test v_none.nonlinear_verbosity isa SciMLLogging.None
        @test v_none.dt_NaN isa SciMLLogging.Silent
        @test v_none.init_NaN isa SciMLLogging.Silent
        @test v_none.alg_switch isa SciMLLogging.Silent
        @test v_none.rosenbrock_no_differential_states isa SciMLLogging.Silent

        # Test Minimal - only critical errors
        @test v_minimal.linear_verbosity isa SciMLLogging.Minimal
        @test v_minimal.nonlinear_verbosity isa SciMLLogging.Minimal
        @test v_minimal.dt_NaN isa SciMLLogging.WarnLevel
        @test v_minimal.init_NaN isa SciMLLogging.WarnLevel
        @test v_minimal.max_iters isa SciMLLogging.WarnLevel
        @test v_minimal.newton_convergence isa SciMLLogging.WarnLevel
        @test v_minimal.near_singular isa SciMLLogging.WarnLevel
        @test v_minimal.alg_switch isa SciMLLogging.Silent
        @test v_minimal.dense_output_saveat isa SciMLLogging.Silent
        @test v_minimal.mismatched_input_output_type isa SciMLLogging.Silent

        # Test Standard - same as default
        @test v_standard.dt_NaN isa SciMLLogging.WarnLevel
        @test v_standard.init_NaN isa SciMLLogging.WarnLevel
        @test v_standard.alg_switch isa SciMLLogging.Silent
        @test v_standard.dense_output_saveat isa SciMLLogging.WarnLevel
        @test v_standard.mismatched_input_output_type isa SciMLLogging.WarnLevel

        # Test Detailed - includes debugging info
        @test v_detailed.linear_verbosity isa SciMLLogging.Detailed
        @test v_detailed.nonlinear_verbosity isa SciMLLogging.Detailed
        @test v_detailed.alg_switch isa SciMLLogging.InfoLevel
        @test v_detailed.dense_output_saveat isa SciMLLogging.InfoLevel
        @test v_detailed.shampine_dt isa SciMLLogging.InfoLevel
        @test v_detailed.jacobian_update isa SciMLLogging.InfoLevel
        @test v_detailed.w_factorization isa SciMLLogging.InfoLevel
        @test v_detailed.convergence_limit isa SciMLLogging.InfoLevel

        # Test All - maximum verbosity
        @test v_all.linear_verbosity isa SciMLLogging.All
        @test v_all.nonlinear_verbosity isa SciMLLogging.All
        @test v_all.alg_switch isa SciMLLogging.InfoLevel
        @test v_all.shampine_dt isa SciMLLogging.InfoLevel
        @test v_all.dense_output_saveat isa SciMLLogging.InfoLevel
        @test v_all.step_rejected isa SciMLLogging.InfoLevel
        @test v_all.step_accepted isa SciMLLogging.InfoLevel
        @test v_all.stiff_detection isa SciMLLogging.InfoLevel
    end

    @testset "Group-level keyword constructors" begin
        v_error = ODEVerbosity(error_control = SciMLLogging.ErrorLevel())
        # Test all error_control fields
        @test v_error.dt_NaN isa SciMLLogging.ErrorLevel
        @test v_error.init_NaN isa SciMLLogging.ErrorLevel
        @test v_error.dense_output_saveat isa SciMLLogging.ErrorLevel
        @test v_error.max_iters isa SciMLLogging.ErrorLevel
        @test v_error.dt_min_unstable isa SciMLLogging.ErrorLevel
        @test v_error.instability isa SciMLLogging.ErrorLevel
        @test v_error.newton_convergence isa SciMLLogging.ErrorLevel
        @test v_error.step_rejected isa SciMLLogging.ErrorLevel
        @test v_error.step_accepted isa SciMLLogging.ErrorLevel
        @test v_error.convergence_limit isa SciMLLogging.ErrorLevel

        v_numerical = ODEVerbosity(numerical = SciMLLogging.Silent())
        # Test all numerical fields
        @test v_numerical.rosenbrock_no_differential_states isa SciMLLogging.Silent
        @test v_numerical.shampine_dt isa SciMLLogging.Silent
        @test v_numerical.unlimited_dt isa SciMLLogging.Silent
        @test v_numerical.dt_epsilon isa SciMLLogging.Silent
        @test v_numerical.stability_check isa SciMLLogging.Silent
        @test v_numerical.near_singular isa SciMLLogging.Silent

        v_performance = ODEVerbosity(performance = SciMLLogging.InfoLevel())
        # Test all performance fields
        @test v_performance.alg_switch isa SciMLLogging.InfoLevel
        @test v_performance.stiff_detection isa SciMLLogging.InfoLevel
        @test v_performance.mismatched_input_output_type isa SciMLLogging.InfoLevel
        @test v_performance.jacobian_update isa SciMLLogging.InfoLevel
        @test v_performance.w_factorization isa SciMLLogging.InfoLevel
        @test v_performance.newton_iterations isa SciMLLogging.InfoLevel
    end

    @testset "Mixed group and individual settings" begin
        v_mixed = ODEVerbosity(
            numerical = SciMLLogging.Silent(),
            shampine_dt = SciMLLogging.WarnLevel(),
            performance = SciMLLogging.InfoLevel()
        )
        # Individual override should take precedence
        @test v_mixed.shampine_dt isa SciMLLogging.WarnLevel
        # Other numerical options should use group setting
        @test v_mixed.rosenbrock_no_differential_states isa SciMLLogging.Silent
        @test v_mixed.unlimited_dt isa SciMLLogging.Silent
        # Performance group setting should apply
        @test v_mixed.alg_switch isa SciMLLogging.InfoLevel
        @test v_mixed.mismatched_input_output_type isa SciMLLogging.InfoLevel
    end

    @testset "Individual keyword arguments" begin
        v_individual = ODEVerbosity(
            dt_NaN = SciMLLogging.ErrorLevel(),
            alg_switch = SciMLLogging.InfoLevel(),
            shampine_dt = SciMLLogging.Silent()
        )
        @test v_individual.dt_NaN isa SciMLLogging.ErrorLevel
        @test v_individual.alg_switch isa SciMLLogging.InfoLevel
        @test v_individual.shampine_dt isa SciMLLogging.Silent
        # Unspecified options should use defaults
        @test v_individual.init_NaN isa SciMLLogging.WarnLevel
        @test v_individual.unlimited_dt isa SciMLLogging.WarnLevel
    end

    @testset "Linear and nonlinear verbosity passthrough" begin
        v_with_solvers = ODEVerbosity(
            linear_verbosity = SciMLLogging.Detailed(),
            nonlinear_verbosity = SciMLLogging.Minimal()
        )
        @test v_with_solvers.linear_verbosity isa SciMLLogging.Detailed
        @test v_with_solvers.nonlinear_verbosity isa SciMLLogging.Minimal

        v_with_solvers2 = ODEVerbosity(
            linear_verbosity = SciMLLogging.None(),
            nonlinear_verbosity = SciMLLogging.All()
        )
        @test v_with_solvers2.linear_verbosity isa SciMLLogging.None
        @test v_with_solvers2.nonlinear_verbosity isa SciMLLogging.All
    end

    @testset "Validation tests" begin
        # Test that invalid group arguments throw errors
        @test_throws ArgumentError ODEVerbosity(error_control = "invalid")
        @test_throws ArgumentError ODEVerbosity(performance = 123)
        @test_throws ArgumentError ODEVerbosity(numerical = :wrong)

        # Test that invalid individual fields throw errors
        @test_throws ArgumentError ODEVerbosity(dt_NaN = "invalid")
        @test_throws ArgumentError ODEVerbosity(unknown_field = SciMLLogging.InfoLevel())
    end

    @testset "Multiple group settings" begin
        v = ODEVerbosity(
            error_control = SciMLLogging.ErrorLevel(),
            performance = SciMLLogging.InfoLevel(),
            numerical = SciMLLogging.Silent()
        )
        @test v.dt_NaN isa SciMLLogging.ErrorLevel
        @test v.alg_switch isa SciMLLogging.InfoLevel
        @test v.shampine_dt isa SciMLLogging.Silent
    end

    @testset "Complex mixed settings" begin
        v = ODEVerbosity(
            error_control = SciMLLogging.WarnLevel(),
            performance = SciMLLogging.InfoLevel(),
            numerical = SciMLLogging.Silent(),
            linear_verbosity = SciMLLogging.Detailed(),
            nonlinear_verbosity = SciMLLogging.Minimal(),
            dt_NaN = SciMLLogging.ErrorLevel(),  # Override specific error_control field
            shampine_dt = SciMLLogging.WarnLevel()  # Override specific numerical field
        )
        # Check overrides took precedence
        @test v.dt_NaN isa SciMLLogging.ErrorLevel
        @test v.shampine_dt isa SciMLLogging.WarnLevel
        # Check other fields follow group settings
        @test v.init_NaN isa SciMLLogging.WarnLevel
        @test v.alg_switch isa SciMLLogging.InfoLevel
        @test v.unlimited_dt isa SciMLLogging.Silent
        # Check solver verbosity
        @test v.linear_verbosity isa SciMLLogging.Detailed
        @test v.nonlinear_verbosity isa SciMLLogging.Minimal
    end

    @testset "Stiff Switching Message" begin
        verb = ODEVerbosity(alg_switch = SciMLLogging.InfoLevel())
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
        tspan = (0.0, 1e-1)
        p = [0.04, 3e7, 1e4]
        prob = ODEProblem(rober, u0, tspan, p)

        @testset "Rosenbrock Solvers" begin
            @testset "Rosenbrock23 with Detailed LinearVerbosity" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.Detailed())
                integrator = init(prob, Rosenbrock23(), verbose = verbose, dt = 1e-3)

                # Check that the cache has a linsolve field
                @test hasproperty(integrator.cache, :linsolve)

                # Check that the linear solver cache has verbose field
                @test hasproperty(integrator.cache.linsolve, :verbose)

                # Verify the verbosity was passed through correctly
                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Detailed())
            end

            @testset "Rosenbrock23 with Minimal LinearVerbosity" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.Minimal())
                integrator = init(prob, Rosenbrock23(), verbose = verbose, dt = 1e-3)

                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Minimal())
            end

            @testset "Rosenbrock23 with None LinearVerbosity" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.None())
                integrator = init(prob, Rosenbrock23(), verbose = verbose, dt = 1e-3)

                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.None())
            end

            @testset "Rosenbrock32 with Detailed LinearVerbosity" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.Detailed())
                integrator = init(prob, Rosenbrock32(), verbose = verbose, dt = 1e-3)

                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Detailed())
            end

            @testset "Rodas4 with All LinearVerbosity" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.All())
                integrator = init(prob, Rodas4(), verbose = verbose, dt = 1e-3)

                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.All())
            end
        end

        @testset "FIRK Solvers (Radau Methods)" begin
            @testset "RadauIIA3 with Detailed LinearVerbosity" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.Detailed())
                integrator = init(prob, RadauIIA3(), verbose = verbose, dt = 1e-3)

                # RadauIIA3 has a linsolve field
                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Detailed())
            end

            @testset "RadauIIA3 with Minimal LinearVerbosity" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.Minimal())
                integrator = init(prob, RadauIIA3(), verbose = verbose, dt = 1e-3)

                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Minimal())
            end

            @testset "RadauIIA5 with Detailed LinearVerbosity (two linear solvers)" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.Detailed())
                integrator = init(prob, RadauIIA5(), verbose = verbose, dt = 1e-3)

                @test integrator.cache.linsolve1.verbose == LinearVerbosity(SciMLLogging.Detailed())
                @test integrator.cache.linsolve2.verbose == LinearVerbosity(SciMLLogging.Detailed())
            end

            @testset "RadauIIA5 with None LinearVerbosity" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.None())
                integrator = init(prob, RadauIIA5(), verbose = verbose, dt = 1e-3)

                @test integrator.cache.linsolve1.verbose == LinearVerbosity(SciMLLogging.None())
                @test integrator.cache.linsolve2.verbose == LinearVerbosity(SciMLLogging.None())
            end

            @testset "RadauIIA9 with All LinearVerbosity (three linear solvers)" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.All())
                integrator = init(prob, RadauIIA9(), verbose = verbose, dt = 1e-3)

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
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.Detailed())
                integrator = init(prob, ImplicitEulerExtrapolation(), verbose = verbose, dt = 1e-3)

                # Check that all linear solvers in the array have the correct verbosity
                for ls in integrator.cache.linsolve
                    @test ls.verbose == LinearVerbosity(SciMLLogging.Detailed())
                end
            end

            @testset "ImplicitEulerExtrapolation with Minimal LinearVerbosity" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.Minimal())
                integrator = init(prob, ImplicitEulerExtrapolation(), verbose = verbose, dt = 1e-3)

                for ls in integrator.cache.linsolve
                    @test ls.verbose == LinearVerbosity(SciMLLogging.Minimal())
                end
            end

            @testset "ImplicitDeuflhardExtrapolation with None LinearVerbosity" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.None())
                integrator = init(prob, ImplicitDeuflhardExtrapolation(), verbose = verbose, dt = 1e-3)

                for ls in integrator.cache.linsolve
                    @test ls.verbose == LinearVerbosity(SciMLLogging.None())
                end
            end

            @testset "ImplicitHairerWannerExtrapolation with All LinearVerbosity" begin
                verbose = ODEVerbosity(linear_verbosity = SciMLLogging.All())
                integrator = init(prob, ImplicitHairerWannerExtrapolation(), verbose = verbose, dt = 1e-3)

                for ls in integrator.cache.linsolve
                    @test ls.verbose == LinearVerbosity(SciMLLogging.All())
                end
            end
        end

        @testset "Preset Verbosity Levels" begin
            @testset "Rosenbrock23 with Standard() preset (default linear_verbosity = Minimal)" begin
                verbose = ODEVerbosity(SciMLLogging.Standard())
                integrator = init(prob, Rosenbrock23(), verbose = verbose, dt = 1e-3)

                # Standard() uses Minimal() for linear_verbosity
                @test integrator.cache.linsolve.verbose == LinearVerbosity(SciMLLogging.Minimal())
            end

            @testset "RadauIIA3 with Detailed() preset" begin
                verbose = ODEVerbosity(SciMLLogging.Detailed())
                integrator = init(prob, RadauIIA3(), verbose = verbose, dt = 1e-3)

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
        tspan = (0.0, 1e-1)
        p = [0.04, 3e7, 1e4]
        prob = ODEProblem(rober, u0, tspan, p)

        @testset "ImplicitEuler with Detailed NonlinearVerbosity" begin
            verbose = ODEVerbosity(nonlinear_verbosity = SciMLLogging.Detailed())
            integrator = init(prob, ImplicitEuler(nlsolve = NonlinearSolveAlg()), verbose = verbose, dt = 1e-3)

            # Check that the cache has an nlsolver field
            @test hasproperty(integrator.cache, :nlsolver)

            # Check that the nlsolver has a cache field
            @test hasproperty(integrator.cache.nlsolver, :cache)

            # Verify the verbosity was passed through correctly
            @test integrator.cache.nlsolver.cache.cache.verbose == NonlinearVerbosity(SciMLLogging.Detailed())
        end

    end
end


