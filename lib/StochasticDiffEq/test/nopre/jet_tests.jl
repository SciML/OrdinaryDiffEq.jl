"""
    JET.jl Static Analysis Tests for StochasticDiffEq.jl

These tests use JET.jl to verify type stability and detect potential
runtime errors in the core SDE solving functionality. The tests focus on
StochasticDiffEq-specific code paths, excluding issues in upstream packages.

For more information on JET.jl, see: https://aviatesk.github.io/JET.jl/stable/
"""

using Test
using StochasticDiffEq
using JET

# Helper to filter JET reports to only show StochasticDiffEq issues
function filter_stochastic_reports(reports)
    return filter(reports) do r
        str_rep = string(r)
        occursin("StochasticDiffEq", str_rep) &&
            !occursin("DiffEqBase", str_rep) &&
            !occursin("SciMLBase", str_rep) &&
            !occursin("OrdinaryDiffEq", str_rep)
    end
end

@testset "JET Static Analysis" begin
    # Define test problems
    function f_iip!(du, u, p, t)
        du[1] = 0.1 * u[1]
        du[2] = 0.2 * u[2]
        du[3] = 0.3 * u[3]
    end
    function g_iip!(du, u, p, t)
        du[1] = 0.2 * u[1]
        du[2] = 0.3 * u[2]
        du[3] = 0.4 * u[3]
    end

    f_scalar(u, p, t) = 0.1u
    g_scalar(u, p, t) = 0.2u

    u0_vec = [1.0, 1.0, 1.0]
    tspan = (0.0, 1.0)

    prob_iip = SDEProblem(f_iip!, g_iip!, u0_vec, tspan)
    prob_scalar = SDEProblem{false}(f_scalar, g_scalar, 1.0, tspan)

    @testset "EM algorithm - in-place" begin
        integrator = init(prob_iip, EM(), dt = 0.01)
        cache = integrator.cache

        # Test perform_step!
        rep = JET.report_call(StochasticDiffEq.perform_step!, (typeof(integrator), typeof(cache)))
        stochastic_issues = filter_stochastic_reports(JET.get_reports(rep))
        @test isempty(stochastic_issues)

        # Test loopheader!
        rep = JET.report_call(StochasticDiffEq.loopheader!, (typeof(integrator),))
        stochastic_issues = filter_stochastic_reports(JET.get_reports(rep))
        @test isempty(stochastic_issues)

        # Test loopfooter!
        rep = JET.report_call(StochasticDiffEq.loopfooter!, (typeof(integrator),))
        stochastic_issues = filter_stochastic_reports(JET.get_reports(rep))
        @test isempty(stochastic_issues)
    end

    @testset "EulerHeun algorithm - in-place" begin
        integrator = init(prob_iip, EulerHeun(), dt = 0.01)
        cache = integrator.cache

        rep = JET.report_call(StochasticDiffEq.perform_step!, (typeof(integrator), typeof(cache)))
        stochastic_issues = filter_stochastic_reports(JET.get_reports(rep))
        @test isempty(stochastic_issues)
    end

    @testset "RKMil algorithm - in-place" begin
        integrator = init(prob_iip, RKMil(), dt = 0.01)
        cache = integrator.cache

        rep = JET.report_call(StochasticDiffEq.perform_step!, (typeof(integrator), typeof(cache)))
        stochastic_issues = filter_stochastic_reports(JET.get_reports(rep))
        @test isempty(stochastic_issues)
    end

    @testset "SOSRI algorithm - in-place (adaptive)" begin
        integrator = init(prob_iip, SOSRI())
        cache = integrator.cache

        rep = JET.report_call(StochasticDiffEq.perform_step!, (typeof(integrator), typeof(cache)))
        stochastic_issues = filter_stochastic_reports(JET.get_reports(rep))
        @test isempty(stochastic_issues)
    end

    @testset "SOSRA algorithm - in-place (adaptive)" begin
        integrator = init(prob_iip, SOSRA())
        cache = integrator.cache

        rep = JET.report_call(StochasticDiffEq.perform_step!, (typeof(integrator), typeof(cache)))
        stochastic_issues = filter_stochastic_reports(JET.get_reports(rep))
        @test isempty(stochastic_issues)
    end

    @testset "SRA1 algorithm - in-place (adaptive)" begin
        integrator = init(prob_iip, SRA1())
        cache = integrator.cache

        rep = JET.report_call(StochasticDiffEq.perform_step!, (typeof(integrator), typeof(cache)))
        stochastic_issues = filter_stochastic_reports(JET.get_reports(rep))
        @test isempty(stochastic_issues)
    end

    @testset "Scalar SDE - out-of-place" begin
        integrator = init(prob_scalar, EM(), dt = 0.01)
        cache = integrator.cache

        rep = JET.report_call(StochasticDiffEq.perform_step!, (typeof(integrator), typeof(cache)))
        stochastic_issues = filter_stochastic_reports(JET.get_reports(rep))
        @test isempty(stochastic_issues)
    end

    @testset "step! function" begin
        integrator = init(prob_iip, EM(), dt = 0.01, save_on = false)

        rep = JET.report_call(step!, (typeof(integrator),))
        stochastic_issues = filter_stochastic_reports(JET.get_reports(rep))
        @test isempty(stochastic_issues)
    end
end
