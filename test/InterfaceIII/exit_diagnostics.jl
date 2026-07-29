using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqBDF
using OrdinaryDiffEqCore: top_error_contributors, error_estimate_residuals
using SciMLBase, StaticArrays, Logging, Test
import OrdinaryDiffEqCore.SciMLLogging as SciMLLogging

# Component 7 is the only stiff one, so it alone drives the explicit method's
# step size down and must be named as the dominant error contributor.
function stiff_component!(du, u, p, t)
    @. du = -0.1 * u
    du[7] = -1.0e5 * u[7]
    return nothing
end
stiff_prob = ODEProblem(stiff_component!, ones(10), (0.0, 10.0))

# Component 4 blows up in finite time while the rest stay bounded.
function blowup_component!(du, u, p, t)
    @. du = u^3
    du[4] = u[4]^5
    return nothing
end
blowup_prob = ODEProblem(blowup_component!, ones(6), (0.0, 10.0))

function solve_logs(args...; kwargs...)
    io = IOBuffer()
    sol = with_logger(SimpleLogger(io)) do
        solve(args...; kwargs...)
    end
    return sol, String(take!(io))
end

@testset "top_error_contributors" begin
    @test top_error_contributors([0.1, -5.0, 3.0, 0.0], 2) == [2, 3]
    @test top_error_contributors([1.0, 2.0], 5) == [2, 1]
    @test top_error_contributors(Float64[], 3) == Int[]
    # NaN outranks every finite residual rather than poisoning the comparisons.
    @test first(top_error_contributors([1.0e9, NaN, 2.0], 3)) == 2
end

@testset "maxiters names the dominant error contributor" begin
    sol, logs = solve_logs(stiff_prob, Tsit5(); maxiters = 500)
    @test sol.retcode == ReturnCode.MaxIters
    @test occursin("Error estimate diagnostics", logs)
    @test occursin("step error estimate EEst", logs)
    @test occursin("atmp[7]", logs)
    @test occursin("u[7]", logs)
    # The exit point is reported so the log alone locates the failure in time.
    @test occursin("Reached t=", logs)
end

@testset "instability exits report the error estimate too" begin
    sol, logs = solve_logs(blowup_prob, Tsit5())
    @test sol.retcode in (ReturnCode.Unstable, ReturnCode.DtLessThanMin)
    @test occursin("atmp[4]", logs)

    sol, logs = solve_logs(blowup_prob, FBDF())
    @test sol.retcode in (ReturnCode.Unstable, ReturnCode.DtLessThanMin)
    @test occursin("atmp[4]", logs)
end

@testset "non-finite residuals are called out" begin
    function nan_component!(du, u, p, t)
        @. du = -u
        du[2] = t > 0.5 ? NaN : -u[2]
        return nothing
    end
    sol, logs = solve_logs(ODEProblem(nan_component!, ones(4), (0.0, 1.0)), Tsit5())
    @test sol.retcode != ReturnCode.Success
    @test occursin("of 4 weighted residuals are non-finite", logs)
    @test occursin("atmp[2] = NaN", logs)
end

@testset "composite and default caches report the active cache" begin
    for alg in (AutoTsit5(Rosenbrock23()),)
        sol, logs = solve_logs(stiff_prob, alg; maxiters = 30)
        @test sol.retcode == ReturnCode.MaxIters
        @test occursin("atmp[7]", logs)
    end
end

@testset "silent verbosity emits nothing" begin
    sol, logs = solve_logs(
        stiff_prob, Tsit5(); maxiters = 500, verbose = SciMLLogging.None()
    )
    @test sol.retcode == ReturnCode.MaxIters
    @test isempty(logs)
end

@testset "degrades gracefully without a stored residual" begin
    # Out-of-place caches build the residuals in a local, so only the scalar
    # estimate is available. The diagnostic must still be emitted.
    oop = ODEProblem(
        (u, p, t) -> SA[-0.1u[1], -0.1u[2], -1.0e5 * u[3]], SA[1.0, 1.0, 1.0],
        (0.0, 10.0)
    )
    sol, logs = solve_logs(oop, Tsit5(); maxiters = 500)
    @test sol.retcode == ReturnCode.MaxIters
    @test occursin("step error estimate EEst", logs)
    @test !occursin("atmp[", logs)

    integ = init(oop, Tsit5())
    @test error_estimate_residuals(integ.cache) === nothing
end

@testset "non-adaptive integration has no error estimate to report" begin
    integ = init(stiff_prob, Tsit5(); adaptive = false, dt = 1.0e-4)
    @test SciMLBase.log_error_estimate(integ) == ""
end
