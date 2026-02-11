using Pkg
Pkg.add("JET")

import OrdinaryDiffEqCore
using OrdinaryDiffEqCore: DEVerbosity
import OrdinaryDiffEqCore.SciMLLogging as SciMLLogging
using JET, Test

@testset "JET Tests" begin
    @test test_package(
        OrdinaryDiffEqCore, target_defined_modules = true, mode = :typo
    ) === nothing broken = true
end

@testset "JET Test DEVerbosity Constructors" begin
    @testset "Default constructor" begin
        @test_opt DEVerbosity()
    end

    @testset "Preset constructors" begin
        @test_opt DEVerbosity(SciMLLogging.None())
        @test_opt DEVerbosity(SciMLLogging.Minimal())
        @test_opt DEVerbosity(SciMLLogging.Standard())
        @test_opt DEVerbosity(SciMLLogging.Detailed())
        @test_opt DEVerbosity(SciMLLogging.All())
    end

    @testset "Group-level keyword constructors" begin
        @test_opt DEVerbosity(error_control = SciMLLogging.ErrorLevel())
        @test_opt DEVerbosity(numerical = SciMLLogging.Silent())
        @test_opt DEVerbosity(performance = SciMLLogging.InfoLevel())
        @test_opt DEVerbosity(error_control = SciMLLogging.WarnLevel())
        @test_opt DEVerbosity(numerical = SciMLLogging.WarnLevel())
        @test_opt DEVerbosity(performance = SciMLLogging.Silent())
    end

    @testset "Individual keyword arguments" begin
        @test_opt DEVerbosity(dt_NaN = SciMLLogging.ErrorLevel())
        @test_opt DEVerbosity(alg_switch = SciMLLogging.InfoLevel())
        @test_opt DEVerbosity(shampine_dt = SciMLLogging.Silent())
        @test_opt DEVerbosity(init_NaN = SciMLLogging.WarnLevel())
    end

    @testset "Linear and nonlinear verbosity" begin
        @test_opt DEVerbosity(linear_verbosity = SciMLLogging.Detailed())
        @test_opt DEVerbosity(nonlinear_verbosity = SciMLLogging.Minimal())
        @test_opt DEVerbosity(linear_verbosity = SciMLLogging.None())
        @test_opt DEVerbosity(nonlinear_verbosity = SciMLLogging.All())
        @test_opt DEVerbosity(
            linear_verbosity = SciMLLogging.Detailed(),
            nonlinear_verbosity = SciMLLogging.Minimal()
        )
        @test_opt DEVerbosity(
            linear_verbosity = SciMLLogging.None(),
            nonlinear_verbosity = SciMLLogging.All()
        )
    end

    @testset "Mixed group and individual settings" begin
        @test_opt DEVerbosity(
            numerical = SciMLLogging.Silent(),
            shampine_dt = SciMLLogging.WarnLevel()
        )
        @test_opt DEVerbosity(
            numerical = SciMLLogging.Silent(),
            shampine_dt = SciMLLogging.WarnLevel(),
            performance = SciMLLogging.InfoLevel()
        )
        @test_opt DEVerbosity(
            error_control = SciMLLogging.WarnLevel(),
            dt_NaN = SciMLLogging.ErrorLevel()
        )
    end

    @testset "Multiple group settings" begin
        @test_opt DEVerbosity(
            error_control = SciMLLogging.ErrorLevel(),
            performance = SciMLLogging.InfoLevel(),
            numerical = SciMLLogging.Silent()
        )
        @test_opt DEVerbosity(
            error_control = SciMLLogging.WarnLevel(),
            performance = SciMLLogging.Silent(),
            numerical = SciMLLogging.WarnLevel()
        )
    end

    @testset "Complex mixed settings" begin
        @test_opt DEVerbosity(
            error_control = SciMLLogging.WarnLevel(),
            performance = SciMLLogging.InfoLevel(),
            numerical = SciMLLogging.Silent(),
            linear_verbosity = SciMLLogging.Detailed(),
            nonlinear_verbosity = SciMLLogging.Minimal(),
            dt_NaN = SciMLLogging.ErrorLevel(),
            shampine_dt = SciMLLogging.WarnLevel()
        )
        @test_opt DEVerbosity(
            error_control = SciMLLogging.ErrorLevel(),
            performance = SciMLLogging.Silent(),
            numerical = SciMLLogging.WarnLevel(),
            linear_verbosity = SciMLLogging.All(),
            nonlinear_verbosity = SciMLLogging.None()
        )
    end
end
