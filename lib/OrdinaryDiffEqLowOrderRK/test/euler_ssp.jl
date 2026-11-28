using OrdinaryDiffEqLowOrderRK
using FastBroadcast: Serial, Threaded
using Polyester
using SciMLBase
using Test

@testset "Euler constructors" begin
    limiter(u, integrator, p, t) = nothing
    F = typeof(limiter)
    for alg in (
            Euler(limiter, limiter), Euler{F, F}(limiter, limiter),
            Euler{F, F}(; stage_limiter! = limiter, step_limiter! = limiter),
        )
        @test alg.stage_limiter! === limiter
        @test alg.step_limiter! === limiter
        @test alg.thread isa Serial
    end
end

@testset "Threaded Euler" begin
    f!(du, u, p, t) = fill!(du, 0.5)
    u0 = fill(1.0, 65536)
    prob = ODEProblem(f!, u0, (0.0, 1.0))
    @test_throws MethodError solve(prob, Euler(thread = :invalid); dt = 0.125)
    serial = solve(prob, Euler(); dt = 0.125)
    integrator = init(prob, Euler(thread = Threaded()); dt = 0.125)
    @test integrator.cache.thread isa Threaded
    threaded = solve!(integrator)
    @test SciMLBase.successful_retcode(threaded)
    @test threaded.u == serial.u
    @test all(u == u0 .+ t * 0.5 for (u, t) in zip(threaded.u, threaded.t))
    @test threaded(0.0625) == u0 .+ 0.0625 * 0.5
    reinit!(integrator, 2 .* u0)
    reused = solve!(integrator)
    @test reused.u[end] == 2 .* u0 .+ 0.5
end

f_ssp = (u, p, t) -> begin
    sin(10t) * u * (1 - u)
end
test_problem_ssp = ODEProblem(f_ssp, 0.1, (0.0, 8.0))
test_problem_ssp_long = ODEProblem(f_ssp, 0.1, (0.0, 1.0e3))

# test SSP coefficient for explicit Euler
alg = Euler()
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqLowOrderRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqLowOrderRK.ssp_coefficient(alg) + 1.0e-3,
    dense = false
)
@test any(sol.u .< 0)
