using OrdinaryDiffEq, OrdinaryDiffEqBDF, OrdinaryDiffEqFIRK, OrdinaryDiffEqRosenbrock
using Test, Logging

# u' = u^2 from u0 = 1 blows up at t = 1, so every solve here ends in retcode Unstable
blowup!(du, u, p, t) = (du[1] = u[1]^2; nothing)
blowup_jac!(J, u, p, t) = (J[1, 1] = 2u[1]; nothing)
blowup(u, p, t) = u .^ 2
blowup_jac(u, p, t) = fill(2u[1], 1, 1)

const U0 = [1.0]
const TSPAN = (0.0, 10.0)

function solve_logs(f, alg; u0 = U0)
    logger = Test.TestLogger(min_level = Logging.Debug)
    sol = with_logger(logger) do
        solve(ODEProblem(f, u0, TSPAN), alg)
    end
    return sol, join((string(r.message) for r in logger.logs), "\n")
end

@testset "Instability diagnostics" begin
    @testset "full diagnostic without an analytic Jacobian" begin
        sol, msg = solve_logs(ODEFunction(blowup!), RadauIIA5())
        @test sol.retcode == ReturnCode.Unstable
        @test occursin("State Analysis", msg)
        @test occursin("Jacobian Analysis", msg)
        @test occursin("Error Analysis", msg)
        @test occursin("J[1,1]", msg)
    end

    # https://github.com/SciML/OrdinaryDiffEq.jl/issues/4155
    @testset "$name analytic Jacobian" for (name, f) in (
            ("in-place", ODEFunction(blowup!; jac = blowup_jac!)),
            ("out-of-place", ODEFunction(blowup; jac = blowup_jac)),
        )
        sol, msg = solve_logs(f, RadauIIA5())
        @test sol.retcode == ReturnCode.Unstable
        @test occursin("Jacobian Analysis", msg)
        @test occursin("J[1,1]", msg)
        @test !occursin("Jacobian analysis unavailable", msg)
    end

    @testset "$(nameof(typeof(alg))) reads its stored Jacobian" for alg in (Rodas5P(), FBDF())
        for f in (ODEFunction(blowup!; jac = blowup_jac!), ODEFunction(blowup; jac = blowup_jac))
            sol, msg = solve_logs(f, alg)
            @test sol.retcode == ReturnCode.Unstable
            @test !occursin("Jacobian analysis unavailable", msg)
        end
    end

    @testset "diagnostic does not inflate njacs" begin
        njac_calls = Ref(0)
        counting_jac!(J, u, p, t) = (njac_calls[] += 1; blowup_jac!(J, u, p, t))
        sol, _ = solve_logs(ODEFunction(blowup!; jac = counting_jac!), RadauIIA5())
        @test sol.retcode == ReturnCode.Unstable
        # the diagnostic evaluates one extra Jacobian that must not reach the stats
        @test njac_calls[] == sol.stats.njacs + 1
    end

    @testset "a failing diagnostic still reports the instability" begin
        njac_calls = Ref(0)
        counting_jac!(J, u, p, t) = (njac_calls[] += 1; blowup_jac!(J, u, p, t))
        solve_logs(ODEFunction(blowup!; jac = counting_jac!), RadauIIA5())
        diagnostic_call = njac_calls[]

        calls = Ref(0)
        function failing_jac!(J, u, p, t)
            calls[] += 1
            calls[] == diagnostic_call && error("deliberate diagnostic failure")
            return blowup_jac!(J, u, p, t)
        end
        sol, msg = solve_logs(ODEFunction(blowup!; jac = failing_jac!), RadauIIA5())
        @test calls[] == diagnostic_call
        @test sol.retcode == ReturnCode.Unstable
        @test occursin("Jacobian analysis unavailable", msg)
        @test occursin("State Analysis", msg)
        @test occursin("Error Analysis", msg)
    end
end
