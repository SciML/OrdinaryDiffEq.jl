using DelayDiffEq
using OrdinaryDiffEqSSPRK
using Test

const LIMITER_CALLS = Ref(0)
limiter!(u, integrator, p, t) = (LIMITER_CALLS[] += 1)

f(du, u, h, p, t) = (du[1] = -0.5 * u[1])
hist(p, t) = [1.0]
prob = DDEProblem(f, [1.0], hist, (0.0, 3.0))

@testset "solve-level step_limiter keyword" begin
    LIMITER_CALLS[] = 0
    sol = solve(prob, MethodOfSteps(SSPRK43()), dt = 0.1; step_limiter = limiter!)
    @test LIMITER_CALLS[] > 0
    @test sol.stats.naccept == LIMITER_CALLS[]
end

@testset "deprecated per-algorithm step_limiter! field is honored" begin
    # `init` depwarns on the deprecated field, which is fatal under
    # `--depwarn=error`, so the deprecated solve runs in a `--depwarn=yes`
    # subprocess; the assertions below then hold in every depwarn mode.
    script = """
    using DelayDiffEq
    using OrdinaryDiffEqSSPRK
    calls = Ref(0)
    limiter!(u, integrator, p, t) = (calls[] += 1)
    f(du, u, h, p, t) = (du[1] = -0.5 * u[1])
    hist(p, t) = [1.0]
    prob = DDEProblem(f, [1.0], hist, (0.0, 3.0))
    sol = solve(prob, MethodOfSteps(SSPRK43(; step_limiter! = limiter!)), dt = 0.1)
    println("LIMITER_RESULT calls=", calls[], " naccept=", sol.stats.naccept)
    """
    out, err = IOBuffer(), IOBuffer()
    proc = run(
        pipeline(
            ignorestatus(`$(Base.julia_cmd()) --depwarn=yes --project=$(dirname(Base.active_project())) -e $script`);
            stdout = out, stderr = err,
        ),
    )
    out_text, err_text = String(take!(out)), String(take!(err))
    m = match(r"LIMITER_RESULT calls=(\d+) naccept=(\d+)", out_text)
    if !success(proc) || m === nothing
        println("subprocess stdout:\n", out_text)
        println("subprocess stderr:\n", err_text)
    end
    @test success(proc)
    @test occursin("`step_limiter!` to the algorithm constructor is deprecated", err_text)
    @test m !== nothing
    if m !== nothing
        calls, naccept = parse.(Int, m.captures)
        @test calls > 0
        @test naccept == calls
    end
end
