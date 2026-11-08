using OrdinaryDiffEqDefault, OrdinaryDiffEqCore, OrdinaryDiffEqVerner, OrdinaryDiffEqBDF, ADTypes, LinearSolve, Test

# A step input makes the problem stiff after t = 0.5, so the composite switches algorithms.
step_rhs!(du, u, p, t) = (du[1] = 1000 * ((t >= 0.5 ? 1.0 : 0.0) - u[1]); nothing)
prob = ODEProblem(step_rhs!, [0.0], (0.0, 1.5))

algs = (
    AutoVern7(FBDF(autodiff = AutoForwardDiff(chunksize = 1), linsolve = LUFactorization())),
    DefaultODEAlgorithm(autodiff = AutoForwardDiff(chunksize = 1)),
)
@testset "reinit! matches a fresh solve: $(nameof(typeof(alg)))" for alg in algs
    fresh = solve(prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
    @test length(unique(fresh.alg_choice)) > 1

    integ = init(prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
    solve!(integ)
    reinit!(integ)
    naccept, nreject = integ.stats.naccept, integ.stats.nreject
    solve!(integ)
    @test integ.sol.t == fresh.t
    @test integ.sol.u == fresh.u
    @test integ.sol.alg_choice == fresh.alg_choice
    @test integ.stats.naccept - naccept == fresh.stats.naccept
    @test integ.stats.nreject - nreject == fresh.stats.nreject
end

import SciMLBase

# Van der Pol switches between the non-stiff and stiff methods several times; the step-size controller's state
# (each method's controller and the error estimate) and the initial dt must be reset by reinit! as at init.
vdp!(du, u, μ, t) = (du[1] = u[2]; du[2] = μ * ((1 - u[1]^2) * u[2] - u[1]); nothing)
vdp = ODEProblem(vdp!, [2.0, 0.0], (0.0, 200.0), 100.0)
@testset "reinit! matches a fresh solve: Van der Pol, DefaultODEAlgorithm" begin
    alg = DefaultODEAlgorithm(autodiff = AutoForwardDiff(chunksize = 1))
    fresh = solve(vdp, alg)
    @test SciMLBase.successful_retcode(fresh)
    @test count(!iszero, diff(fresh.alg_choice)) > 1     # switches back and forth
    integ = init(vdp, alg)
    solve!(integ)
    reinit!(integ)
    @test integ.controller_cache.EEst == 1
    naccept = integ.stats.naccept
    solve!(integ)
    @test SciMLBase.successful_retcode(integ.sol)
    @test integ.sol.t == fresh.t
    @test integ.sol.u == fresh.u
    @test integ.sol.alg_choice == fresh.alg_choice
    @test integ.stats.naccept - naccept == fresh.stats.naccept
end

# reinit! keeps a dt that a callback's initialization sets, and does not estimate one after a callback's
# initialization terminated the integration
decay!(du, u, p, t) = (du[1] = -u[1]; nothing)
decay = ODEProblem(decay!, [1.0], (0.0, 1.0))
@testset "reinit! with callback initializations that set dt or terminate" begin
    inits = Ref(0)
    set_dt = SciMLBase.DiscreteCallback(
        (u, t, integ) -> false, integ -> nothing;
        initialize = (c, u, t, integ) -> (inits[] += 1; inits[] > 1 && (integ.dt = 0.0123); nothing)
    )
    integ = init(decay, DefaultODEAlgorithm(); callback = set_dt)
    solve!(integ)
    reinit!(integ)
    @test integ.dt == 0.0123
    inits[] = 0
    stop = SciMLBase.DiscreteCallback(
        (u, t, integ) -> false, integ -> nothing;
        initialize = (c, u, t, integ) -> (inits[] += 1; inits[] > 1 && SciMLBase.terminate!(integ); nothing)
    )
    integ = init(decay, DefaultODEAlgorithm(); callback = stop)
    solve!(integ)
    reinit!(integ)
    @test integ.sol.retcode == SciMLBase.ReturnCode.Terminated
    solve!(integ)
    @test integ.sol.retcode == SciMLBase.ReturnCode.Terminated
end

# u = 0 before the step and 1 - exp(-1000 (t - 0.5)) after it
step_exact(t) = t <= 0.5 ? 0.0 : -expm1(-1000 * (t - 0.5))
@testset "reinit! and solve! against the exact solution: $(nameof(typeof(alg)))" for alg in algs
    ts = [0.0, 0.25, 0.5, 0.5001, 0.501, 1.5]
    integ = init(prob, alg; abstol = 1.0e-10, reltol = 1.0e-10, tstops = [0.5])
    fresh = solve(prob, alg; abstol = 1.0e-10, reltol = 1.0e-10, tstops = [0.5])
    for run in 1:3
        run > 1 && reinit!(integ)
        solve!(integ)
        @test SciMLBase.successful_retcode(integ.sol)
        @test all(t -> isapprox(integ.sol(t)[1], step_exact(t); atol = 1.0e-7, rtol = 1.0e-6), ts)
        @test integ.sol.t == fresh.t
    end
    # a new time span, as a fresh solve of the remade problem
    reinit!(integ; t0 = 0.25, tf = 1.0)
    solve!(integ)
    fresh = solve(remake(prob; tspan = (0.25, 1.0)), alg; abstol = 1.0e-10, reltol = 1.0e-10, tstops = [0.5])
    @test integ.sol.t == fresh.t
    @test integ.sol.u == fresh.u
    @test all(t -> isapprox(integ.sol(t)[1], step_exact(t); atol = 1.0e-7, rtol = 1.0e-6), (0.25, 0.5, 0.501, 1.0))
end

member_caches(cache::OrdinaryDiffEqCore.CompositeCache) = collect(pairs(cache.caches))
function member_caches(cache::OrdinaryDiffEqCore.DefaultCache)
    members = Symbol.(:cache, 1:6)
    return [k => getfield(cache, m) for (k, m) in enumerate(members) if isdefined(cache, m)]
end
function nlsolver_states(cache)
    return [
        k => (sub.nlsolver.prev_θ, sub.nlsolver.cache.firstcall)
            for (k, sub) in member_caches(cache) if hasfield(typeof(sub), :nlsolver)
    ]
end

# With more than 50 states, the default algorithm's stiff member is FBDF.
step_rhs_60!(du, u, p, t) = (du .= 1000 .* ((t >= 0.5 ? 1.0 : 0.0) .- u); nothing)
prob_60 = ODEProblem(step_rhs_60!, zeros(60), (0.0, 1.5))
state_cases = (
    "AutoVern7(FBDF())" => (prob, algs[1]),
    "DefaultODEAlgorithm, 60 states" => (prob_60, algs[2]),
)
@testset "state after reinit! equals the state after init: $name" for (name, (prob, alg)) in state_cases
    fresh_sol = solve(prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
    integ = init(prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
    solve!(integ)
    @test !isempty(nlsolver_states(integ.cache))
    @test any(((k, s),) -> s != (1.0, true), nlsolver_states(integ.cache))
    reinit!(integ)
    fresh = init(prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
    if integ.cache isa OrdinaryDiffEqCore.DefaultCache
        for (k, _) in member_caches(integ.cache)
            OrdinaryDiffEqCore.init_ith_default_cache(fresh.cache, fresh.alg.algs, k)
        end
    end
    @test integ.dt == fresh.dt
    @test integ.cache.current == fresh.cache.current
    @test integ.eigen_est == fresh.eigen_est
    @test integ.controller_cache.EEst == fresh.controller_cache.EEst
    AS, AS0 = integ.alg.choice_function, fresh.alg.choice_function
    @test (AS.count, AS.successive_switches, AS.is_stiffalg, AS.current) ==
        (AS0.count, AS0.successive_switches, AS0.is_stiffalg, AS0.current)
    @test nlsolver_states(integ.cache) == nlsolver_states(fresh.cache)
    solve!(integ)
    @test integ.sol.t == fresh_sol.t
    @test integ.sol.u == fresh_sol.u
    @test integ.sol.alg_choice == fresh_sol.alg_choice
    @test all(t -> isapprox(integ.sol(t)[end], step_exact(t); atol = 1.0e-6, rtol = 1.0e-5), (0.25, 0.501, 1.5))
end

@testset "reinit! of one integrator leaves another one built from the same algorithm alone" begin
    alg = DefaultODEAlgorithm(autodiff = AutoForwardDiff(chunksize = 1))
    a = init(vdp, alg)
    b = init(vdp, alg)
    solve!(a)
    solve!(b)
    AS = b.alg.choice_function
    before = (AS.count, AS.successive_switches, AS.is_stiffalg, AS.current, b.cache.current)
    reinit!(a)
    @test (AS.count, AS.successive_switches, AS.is_stiffalg, AS.current, b.cache.current) == before
end

@testset "reinit! keeps the switching and nonlinear solver state with reinit_cache = false" begin
    integ = init(prob_60, algs[2]; abstol = 1.0e-8, reltol = 1.0e-8)
    solve!(integ)
    AS = integ.alg.choice_function
    before = (AS.count, AS.successive_switches, AS.is_stiffalg, AS.current)
    nls = nlsolver_states(integ.cache)
    @test !isempty(nls)
    reinit!(integ; reinit_cache = false)
    @test (AS.count, AS.successive_switches, AS.is_stiffalg, AS.current) == before
    @test nlsolver_states(integ.cache) == nls
end

# Robertson in mass-matrix form; u1 + u2 + u3 = 1 holds exactly for the true solution.
function rober!(du, u, p, t)
    y1, y2, y3 = u
    du[1] = -0.04y1 + 1.0e4 * y2 * y3
    du[2] = 0.04y1 - 3.0e7 * y2^2 - 1.0e4 * y2 * y3
    du[3] = y1 + y2 + y3 - 1
    return nothing
end
rober = ODEProblem(
    ODEFunction(rober!; mass_matrix = [1.0 0 0; 0 1 0; 0 0 0]), [1.0, 0.0, 0.0], (0.0, 1.0e5)
)
@testset "reinit! matches a fresh solve: mass-matrix DAE, DefaultODEAlgorithm" begin
    alg = DefaultODEAlgorithm(autodiff = AutoForwardDiff(chunksize = 1))
    fresh = solve(rober, alg; abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(fresh)
    integ = init(rober, alg; abstol = 1.0e-8, reltol = 1.0e-8)
    solve!(integ)
    reinit!(integ)
    solve!(integ)
    @test SciMLBase.successful_retcode(integ.sol)
    @test integ.sol.t == fresh.t
    @test integ.sol.u == fresh.u
    @test all(u -> isapprox(sum(u), 1; atol = 1.0e-8), integ.sol.u)
end
