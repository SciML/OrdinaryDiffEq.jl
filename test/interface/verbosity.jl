using OrdinaryDiffEqCore
using OrdinaryDiffEqCore: ODEVerbosity, option_group, group_options
using ODEProblemLibrary: prob_ode_vanderpol_stiff
using SciMLLogging
using Test

@testset "ODEVerbosity Tests" begin
    @testset "Default constructor" begin
        v1 = ODEVerbosity()
        @test v1 isa ODEVerbosity
        @test v1.dt_NaN isa SciMLLogging.WarnLevel
        @test v1.init_NaN isa SciMLLogging.WarnLevel
        @test v1.dense_output_saveat isa SciMLLogging.WarnLevel
        @test v1.alg_switch isa SciMLLogging.WarnLevel
        @test v1.mismatched_input_output_type isa SciMLLogging.WarnLevel
        @test v1.rosenbrock_no_differential_states isa SciMLLogging.WarnLevel
        @test v1.shampine_dt isa SciMLLogging.WarnLevel
        @test v1.unlimited_dt isa SciMLLogging.WarnLevel
        @test v1.linear_verbosity isa SciMLLogging.AbstractVerbosityPreset
        @test v1.nonlinear_verbosity isa SciMLLogging.AbstractVerbosityPreset
    end

    @testset "ODEVerbosity preset constructors" begin
        v_none = ODEVerbosity(SciMLLogging.None())
        v_all = ODEVerbosity(SciMLLogging.All())
        v_minimal = ODEVerbosity(SciMLLogging.Minimal())
        v_standard = ODEVerbosity(SciMLLogging.Standard())
        v_detailed = ODEVerbosity(SciMLLogging.Detailed())

        @test v_none.dt_NaN isa SciMLLogging.Silent
        @test v_none.init_NaN isa SciMLLogging.Silent
        @test v_none.alg_switch isa SciMLLogging.Silent
        @test v_none.rosenbrock_no_differential_states isa SciMLLogging.Silent

        @test v_minimal.dt_NaN isa SciMLLogging.ErrorLevel
        @test v_minimal.init_NaN isa SciMLLogging.ErrorLevel
        @test v_minimal.alg_switch isa SciMLLogging.Silent
        @test v_minimal.dense_output_saveat isa SciMLLogging.Silent

        @test v_standard.dt_NaN isa SciMLLogging.WarnLevel
        @test v_standard.init_NaN isa SciMLLogging.WarnLevel
        @test v_standard.alg_switch isa SciMLLogging.WarnLevel

        @test v_detailed.alg_switch isa SciMLLogging.InfoLevel
        @test v_detailed.dense_output_saveat isa SciMLLogging.InfoLevel
        @test v_detailed.shampine_dt isa SciMLLogging.InfoLevel

        @test v_all.alg_switch isa SciMLLogging.InfoLevel
        @test v_all.shampine_dt isa SciMLLogging.InfoLevel
        @test v_all.dense_output_saveat isa SciMLLogging.InfoLevel
    end

    @testset "Group-level keyword constructors" begin
        v_error = ODEVerbosity(error_control = ErrorLevel())
        @test v_error.dt_NaN isa SciMLLogging.ErrorLevel
        @test v_error.init_NaN isa SciMLLogging.ErrorLevel
        @test v_error.dense_output_saveat isa SciMLLogging.ErrorLevel

        v_numerical = ODEVerbosity(numerical = Silent())
        @test v_numerical.rosenbrock_no_differential_states isa SciMLLogging.Silent
        @test v_numerical.shampine_dt isa SciMLLogging.Silent
        @test v_numerical.unlimited_dt isa SciMLLogging.Silent

        v_performance = ODEVerbosity(performance = InfoLevel())
        @test v_performance.alg_switch isa SciMLLogging.InfoLevel
        @test v_performance.mismatched_input_output_type isa SciMLLogging.InfoLevel
    end

    @testset "Mixed group and individual settings" begin
        v_mixed = ODEVerbosity(
            numerical = Silent(),
            shampine_dt = WarnLevel(),
            performance = InfoLevel()
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
            dt_NaN = ErrorLevel(),
            alg_switch = InfoLevel(),
            shampine_dt = Silent()
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

    @testset "Group classification functions" begin
        @test option_group(:dt_NaN) == :error_control
        @test option_group(:init_NaN) == :error_control
        @test option_group(:dense_output_saveat) == :error_control
        @test option_group(:alg_switch) == :performance
        @test option_group(:mismatched_input_output_type) == :performance
        @test option_group(:rosenbrock_no_differential_states) == :numerical
        @test option_group(:shampine_dt) == :numerical
        @test option_group(:unlimited_dt) == :numerical

        # Test error for unknown option
        @test_throws ErrorException option_group(:unknown_option)
    end

    @testset "Group options function" begin
        v = ODEVerbosity(numerical = WarnLevel())
        numerical_opts = group_options(v, :numerical)
        @test numerical_opts isa NamedTuple
        @test :rosenbrock_no_differential_states in keys(numerical_opts)
        @test :shampine_dt in keys(numerical_opts)
        @test :unlimited_dt in keys(numerical_opts)
        @test numerical_opts.rosenbrock_no_differential_states isa SciMLLogging.WarnLevel
        @test numerical_opts.shampine_dt isa SciMLLogging.WarnLevel
        @test numerical_opts.unlimited_dt isa SciMLLogging.WarnLevel

        error_opts = group_options(v, :error_control)
        @test :dt_NaN in keys(error_opts)
        @test :init_NaN in keys(error_opts)
        @test :dense_output_saveat in keys(error_opts)

        performance_opts = group_options(v, :performance)
        @test :alg_switch in keys(performance_opts)
        @test :mismatched_input_output_type in keys(performance_opts)

        # Test error for unknown group
        @test_throws ErrorException group_options(v, :unknown_group)
    end

    @testset "Group getproperty access" begin
        v = ODEVerbosity()

        # Test getting groups returns NamedTuples
        error_group = v.error_control
        performance_group = v.performance
        numerical_group = v.numerical

        @test error_group isa NamedTuple
        @test performance_group isa NamedTuple
        @test numerical_group isa NamedTuple

        # Test correct keys are present
        @test :dt_NaN in keys(error_group)
        @test :init_NaN in keys(error_group)
        @test :dense_output_saveat in keys(error_group)

        @test :alg_switch in keys(performance_group)
        @test :mismatched_input_output_type in keys(performance_group)

        @test :rosenbrock_no_differential_states in keys(numerical_group)
        @test :shampine_dt in keys(numerical_group)
        @test :unlimited_dt in keys(numerical_group)

        # Test values are AbstractMessageLevel types
        @test error_group.dt_NaN isa SciMLLogging.AbstractMessageLevel
        @test performance_group.alg_switch isa SciMLLogging.AbstractMessageLevel
        @test numerical_group.shampine_dt isa SciMLLogging.AbstractMessageLevel

        # Individual field access should still work
        @test v.dt_NaN isa SciMLLogging.WarnLevel
        @test v.alg_switch isa SciMLLogging.WarnLevel
        @test v.shampine_dt isa SciMLLogging.WarnLevel
    end

    @testset "Argument validation" begin
        # Test invalid error_control type
        @test_throws ArgumentError ODEVerbosity(error_control = "invalid")
        @test_throws ArgumentError ODEVerbosity(performance = 123)
        @test_throws ArgumentError ODEVerbosity(numerical = :symbol)

        # Test unknown keyword argument
        @test_throws ArgumentError ODEVerbosity(unknown_field = WarnLevel())

        # Test invalid value for individual field
        @test_throws ArgumentError ODEVerbosity(dt_NaN = "not_a_level")
    end

    @testset "All error control fields" begin
        v = ODEVerbosity(error_control = InfoLevel())
        @test v.dt_NaN isa SciMLLogging.InfoLevel
        @test v.init_NaN isa SciMLLogging.InfoLevel
        @test v.dense_output_saveat isa SciMLLogging.InfoLevel
    end

    @testset "All performance fields" begin
        v = ODEVerbosity(performance = ErrorLevel())
        @test v.alg_switch isa SciMLLogging.ErrorLevel
        @test v.mismatched_input_output_type isa SciMLLogging.ErrorLevel
    end

    @testset "All numerical fields" begin
        v = ODEVerbosity(numerical = InfoLevel())
        @test v.rosenbrock_no_differential_states isa SciMLLogging.InfoLevel
        @test v.shampine_dt isa SciMLLogging.InfoLevel
        @test v.unlimited_dt isa SciMLLogging.InfoLevel
    end

    @testset "Multiple group settings" begin
        v = ODEVerbosity(
            error_control = ErrorLevel(),
            performance = InfoLevel(),
            numerical = Silent()
        )
        @test v.dt_NaN isa SciMLLogging.ErrorLevel
        @test v.alg_switch isa SciMLLogging.InfoLevel
        @test v.shampine_dt isa SciMLLogging.Silent
    end

    @testset "Complex mixed settings" begin
        v = ODEVerbosity(
            error_control = WarnLevel(),
            performance = InfoLevel(),
            numerical = Silent(),
            linear_verbosity = SciMLLogging.Detailed(),
            nonlinear_verbosity = SciMLLogging.Minimal(),
            dt_NaN = ErrorLevel(),  # Override specific error_control field
            shampine_dt = WarnLevel()  # Override specific numerical field
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
        verb = ODEVerbosity(performance = ODEPerformanceVerbosity(alg_switch = Verbosity.Info()))
        solve(prob_ode_vanderpol_stiff, AutoTsit5(Rodas5()), verbose = verb)
    end 
end


