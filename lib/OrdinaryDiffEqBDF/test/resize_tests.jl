using OrdinaryDiffEqBDF, LinearAlgebra, OffsetArrays, Test
using OrdinaryDiffEqNonlinearSolve: NLFunctional

const A = [-1.0 0.35 0.1; 0.2 -1.4 0.25; 0.15 0.3 -1.2]
const U0 = [1.0, -0.5]
const UNEW = 0.3
const TG = 0.5
const TS = 1.0
const TF = 2.0

f!(du, u, p, t) = (mul!(du, view(A, 1:length(u), 1:length(u)), u); nothing)

function exact(t)
    u = exp(A[1:2, 1:2] * min(t, TG)) * U0
    t <= TG && return u
    u3 = exp(A * (min(t, TS) - TG)) * vcat(u, UNEW)
    t <= TS && return u3
    return exp(A[1:2, 1:2] * (t - TS)) * u3[1:2]
end

function resize_callback!(disc, restart_orders)
    grew = Ref(false)
    shrank = Ref(false)
    function condition(u, t, integrator)
        return (!grew[] && t >= TG) || (!shrank[] && t >= TS)
    end
    function affect(integrator)
        n = grew[] ? 2 : 3
        resize!(integrator, n)
        n == 3 && (integrator.u[3] = UNEW)
        push!(restart_orders, integrator.cache.order)
        grew[] ? (shrank[] = true) : (grew[] = true)
        return derivative_discontinuity!(integrator, disc)
    end
    return DiscreteCallback(condition, affect; save_positions = (false, false))
end

@testset "BDF state resizing" begin
    prob = ODEProblem(f!, copy(U0), (0.0, TF))

    for (alg, restart_order) in ((FBDF(), 1), (QNDF(), 1), (MOOSE234(), 2)), disc in (true, false)
        restart_orders = Int[]
        sol = solve(
            prob, alg; callback = resize_callback!(disc, restart_orders),
            tstops = [TG, TS], reltol = 1.0e-10, abstol = 1.0e-12
        )
        @test sol(0.75) ≈ exact(0.75) rtol = 1.0e-7 atol = 1.0e-10
        @test sol.u[end] ≈ exact(TF) rtol = 1.0e-7 atol = 1.0e-10
        @test restart_orders == fill(restart_order, 2)
    end
end

# Shrinking a vector with `resize!` keeps its trailing memory, so growing it back exposes the old
# values. QNDF keeps its difference history across a resize when the derivative discontinuity is
# cleared, so the rows a resize adds to that history must not carry such values.
@testset "QNDF history rows added by a resize" begin
    prob = ODEProblem(f!, [U0; UNEW], (0.0, TF))
    integrator = init(prob, QNDF(); reltol = 1.0e-10, abstol = 1.0e-12)
    for _ in 1:10
        step!(integrator)
    end
    history = (integrator.cache.D, integrator.cache.Dtmp, integrator.cache.prevD)
    foreach(d -> foreach(x -> x[3] = NaN, d), history)
    u3 = integrator.u[3]
    resize!(integrator, 2)
    resize!(integrator, 3)
    integrator.u[3] = u3
    derivative_discontinuity!(integrator, false)
    @test all(x[3] == 0 for d in history for x in d)
    solve!(integrator)
    @test integrator.sol.retcode == ReturnCode.Success
    @test integrator.u ≈ exp(A * TF) * [U0; UNEW] rtol = 1.0e-7
end

decay!(du, u, p, t) = (du .= -u; nothing)

@testset "QNDF history resize with offset axes ($offset)" for offset in (-2, 2)
    prob = ODEProblem(decay!, OffsetArray(collect(1.0:5.0), offset), (0.0, 1.0))
    integrator = init(prob, QNDF(nlsolve = NLFunctional()); reltol = 1.0e-8, abstol = 1.0e-10)
    for _ in 1:5
        step!(integrator)
    end
    cache = integrator.cache
    history = [x for d in (cache.D, cache.Dtmp, cache.prevD) for x in d]
    kept = map(x -> copy(parent(x)), history)
    resize!(integrator, 6)
    @test all(parent(x)[1:5] == k for (x, k) in zip(history, kept))
    @test all(x[end] == 0 for x in history)
    integrator.u[end] = 0.0
    derivative_discontinuity!(integrator, false)
    solve!(integrator)
    @test integrator.sol.retcode == ReturnCode.Success
    @test parent(integrator.u) ≈ [exp(-1.0) .* (1:5); 0.0] rtol = 1.0e-6
end
