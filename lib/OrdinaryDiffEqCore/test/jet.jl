using Pkg
Pkg.add("JET")

import OrdinaryDiffEqCore
using OrdinaryDiffEqCore: ODEVerbosity
import OrdinaryDiffEqCore.SciMLLogging as SciMLLogging
using JET, Test

@testset "JET Tests" begin
    @test test_package(
        OrdinaryDiffEqCore, target_defined_modules = true, mode = :typo
    ) === nothing broken = true
end

@testset "JET Test ODEVerbosity Constructors" begin
    @testset "Default constructor" begin
        @test_opt ODEVerbosity()
    end

    @testset "Preset constructors" begin
        @test_opt ODEVerbosity(SciMLLogging.None())
        @test_opt ODEVerbosity(SciMLLogging.Minimal())
        @test_opt ODEVerbosity(SciMLLogging.Standard())
        @test_opt ODEVerbosity(SciMLLogging.Detailed())
        @test_opt ODEVerbosity(SciMLLogging.All())
    end

    @testset "Group-level keyword constructors" begin
        @test_opt ODEVerbosity(error_control = SciMLLogging.ErrorLevel())
        @test_opt ODEVerbosity(numerical = SciMLLogging.Silent())
        @test_opt ODEVerbosity(performance = SciMLLogging.InfoLevel())
        @test_opt ODEVerbosity(error_control = SciMLLogging.WarnLevel())
        @test_opt ODEVerbosity(numerical = SciMLLogging.WarnLevel())
        @test_opt ODEVerbosity(performance = SciMLLogging.Silent())
    end

    @testset "Individual keyword arguments" begin
        @test_opt ODEVerbosity(dt_NaN = SciMLLogging.ErrorLevel())
        @test_opt ODEVerbosity(alg_switch = SciMLLogging.InfoLevel())
        @test_opt ODEVerbosity(shampine_dt = SciMLLogging.Silent())
        @test_opt ODEVerbosity(init_NaN = SciMLLogging.WarnLevel())
    end

    @testset "Linear and nonlinear verbosity" begin
        @test_opt ODEVerbosity(linear_verbosity = SciMLLogging.Detailed())
        @test_opt ODEVerbosity(nonlinear_verbosity = SciMLLogging.Minimal())
        @test_opt ODEVerbosity(linear_verbosity = SciMLLogging.None())
        @test_opt ODEVerbosity(nonlinear_verbosity = SciMLLogging.All())
        @test_opt ODEVerbosity(
            linear_verbosity = SciMLLogging.Detailed(),
            nonlinear_verbosity = SciMLLogging.Minimal()
        )
        @test_opt ODEVerbosity(
            linear_verbosity = SciMLLogging.None(),
            nonlinear_verbosity = SciMLLogging.All()
        )
    end

    @testset "Mixed group and individual settings" begin
        @test_opt ODEVerbosity(
            numerical = SciMLLogging.Silent(),
            shampine_dt = SciMLLogging.WarnLevel()
        )
        @test_opt ODEVerbosity(
            numerical = SciMLLogging.Silent(),
            shampine_dt = SciMLLogging.WarnLevel(),
            performance = SciMLLogging.InfoLevel()
        )
        @test_opt ODEVerbosity(
            error_control = SciMLLogging.WarnLevel(),
            dt_NaN = SciMLLogging.ErrorLevel()
        )
    end

    @testset "Multiple group settings" begin
        @test_opt ODEVerbosity(
            error_control = SciMLLogging.ErrorLevel(),
            performance = SciMLLogging.InfoLevel(),
            numerical = SciMLLogging.Silent()
        )
        @test_opt ODEVerbosity(
            error_control = SciMLLogging.WarnLevel(),
            performance = SciMLLogging.Silent(),
            numerical = SciMLLogging.WarnLevel()
        )
    end

    @testset "Complex mixed settings" begin
        @test_opt ODEVerbosity(
            error_control = SciMLLogging.WarnLevel(),
            performance = SciMLLogging.InfoLevel(),
            numerical = SciMLLogging.Silent(),
            linear_verbosity = SciMLLogging.Detailed(),
            nonlinear_verbosity = SciMLLogging.Minimal(),
            dt_NaN = SciMLLogging.ErrorLevel(),
            shampine_dt = SciMLLogging.WarnLevel()
        )
        @test_opt ODEVerbosity(
            error_control = SciMLLogging.ErrorLevel(),
            performance = SciMLLogging.Silent(),
            numerical = SciMLLogging.WarnLevel(),
            linear_verbosity = SciMLLogging.All(),
            nonlinear_verbosity = SciMLLogging.None()
        )
    end
end
