using DelayDiffEq, Test
using OrdinaryDiffEqCore: DEVerbosity
using SciMLLogging: None, Minimal, Standard, Detailed, All
using SciMLBase: ReturnCode

# Simple delay problem for testing
function f_dde(du, u, h, p, t)
    du[1] = -u[1] + h(p, t - 1.0)[1]
    return nothing
end

h(p, t) = [1.0]
u0 = [1.0]
tspan = (0.0, 2.0)
prob = DDEProblem(f_dde, u0, h, tspan; constant_lags = [1.0])

@testset "Verbosity Tests" begin
    # Test default verbosity (should work)
    @testset "Default verbosity" begin
        sol = solve(prob; saveat = 0.1)
        @test sol.retcode == ReturnCode.Success
    end

    # Test that Bool verbose is no longer supported (v7 breaking change)
    @testset "verbose=true throws" begin
        @test_throws ArgumentError solve(prob; verbose = true, saveat = 0.1)
    end

    @testset "verbose=false throws" begin
        @test_throws ArgumentError solve(prob; verbose = false, saveat = 0.1)
    end

    # Test with DEVerbosity
    @testset "DEVerbosity()" begin
        sol = solve(prob; verbose = DEVerbosity(), saveat = 0.1)
        @test sol.retcode == ReturnCode.Success
    end

    # Test with preset verbosity levels
    @testset "Preset verbosity levels" begin
        sol = solve(prob; verbose = DEVerbosity(None()), saveat = 0.1)
        @test sol.retcode == ReturnCode.Success

        sol = solve(prob; verbose = DEVerbosity(Minimal()), saveat = 0.1)
        @test sol.retcode == ReturnCode.Success

        sol = solve(prob; verbose = DEVerbosity(Standard()), saveat = 0.1)
        @test sol.retcode == ReturnCode.Success

        sol = solve(prob; verbose = DEVerbosity(Detailed()), saveat = 0.1)
        @test sol.retcode == ReturnCode.Success

        sol = solve(prob; verbose = DEVerbosity(All()), saveat = 0.1)
        @test sol.retcode == ReturnCode.Success
    end
end
