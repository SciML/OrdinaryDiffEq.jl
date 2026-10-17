using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using DiffEqDevTools
using DiffEqBase
using LineSearches
using Test

using ODEProblemLibrary: prob_ode_lorenz, prob_ode_orego

for prob in (prob_ode_lorenz, prob_ode_orego)
    sol1 = solve(prob, Trapezoid(), reltol = 1.0e-12, abstol = 1.0e-12)
    @test sol1.retcode == SciMLBase.ReturnCode.Success
    sol2 = solve(
        prob, Trapezoid(nlsolve = NLNewton(relax = BackTracking())),
        reltol = 1.0e-12, abstol = 1.0e-12
    )
    @test sol2.retcode == SciMLBase.ReturnCode.Success
    @test sol2.stats.nf <= sol1.stats.nf + 20
end

# Checked here because older LineSearches releases accept a positive initial slope silently.
struct SlopeRecorder{L}
    linesearch::L
    slopes::Vector{Float64}
end
function (r::SlopeRecorder)(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
    push!(r.slopes, dϕ0)
    return r.linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
end

for prob in (prob_ode_lorenz, prob_ode_orego)
    recorder = SlopeRecorder(BackTracking(), Float64[])
    sol = solve(
        prob, Trapezoid(nlsolve = NLNewton(relax = recorder)),
        reltol = 1.0e-12, abstol = 1.0e-12
    )
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test !isempty(recorder.slopes)
    @test count(>(0), recorder.slopes) == 0
end

@testset "In-place line-search accuracy" begin
    function riccati!(du, u, p, t)
        du[1] = -10 * u[1]^2
    end
    prob = ODEProblem(riccati!, [1.0], (0.0, 1.0))
    recorder = SlopeRecorder(BackTracking(), Float64[])
    sol = solve(
        prob, Trapezoid(nlsolve = NLNewton(relax = recorder)),
        reltol = 1.0e-12, abstol = 1.0e-14, saveat = 0.0:0.1:1.0
    )
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test !isempty(recorder.slopes)
    for (t, u) in zip(sol.t, sol.u)
        @test isapprox(u[1], inv(1 + 10t); rtol = 1.0e-7, atol = 1.0e-9)
    end
end
