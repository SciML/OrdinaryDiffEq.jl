using OrdinaryDiffEqDefault, OrdinaryDiffEqVerner, OrdinaryDiffEqBDF, ADTypes, LinearSolve, Test

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

# Van der Pol switches between the non-stiff and stiff methods several times; the step-size controller's state
# (each method's controller and the error estimate) must be reset by reinit! as at init.
vdp!(du, u, μ, t) = (du[1] = u[2]; du[2] = μ * ((1 - u[1]^2) * u[2] - u[1]); nothing)
vdp = ODEProblem(vdp!, [2.0, 0.0], (0.0, 200.0), 100.0)
@testset "reinit! matches a fresh solve: Van der Pol, DefaultODEAlgorithm" begin
    alg = DefaultODEAlgorithm(autodiff = AutoForwardDiff(chunksize = 1))
    fresh = solve(vdp, alg; save_everystep = false)
    integ = init(vdp, alg; save_everystep = false)
    solve!(integ)
    reinit!(integ)
    @test integ.controller_cache.EEst == 1
    naccept = integ.stats.naccept
    solve!(integ)
    @test integ.sol.t == fresh.t
    @test integ.sol.u == fresh.u
    @test integ.stats.naccept - naccept == fresh.stats.naccept
end

