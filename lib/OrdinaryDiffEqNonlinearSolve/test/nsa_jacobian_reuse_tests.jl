using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg, NLNewton
using NonlinearSolve: NewtonRaphson, TrustRegion
using ADTypes, LinearAlgebra, SciMLBase
using Test

function rober!(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    du[3] = k₂ * y₂^2
    return nothing
end
const JAC_CALLS = Ref(0)
function rober_jac!(J, u, p, t)
    JAC_CALLS[] += 1
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    J[1, 1] = -k₁
    J[1, 2] = k₃ * y₃
    J[1, 3] = k₃ * y₂
    J[2, 1] = k₁
    J[2, 2] = -2k₂ * y₂ - k₃ * y₃
    J[2, 3] = -k₃ * y₂
    J[3, 1] = 0
    J[3, 2] = 2k₂ * y₂
    J[3, 3] = 0
    return nothing
end
f = ODEFunction(rober!; jac = rober_jac!)
prob = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1.0e5), [0.04, 3.0e7, 1.0e4])
nsa = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff()))
refsol = solve(prob, FBDF(); reltol = 1.0e-12, abstol = 1.0e-14)

@testset "dt-only W updates reuse the stored Jacobian" begin
    for alg in (FBDF(nlsolve = nsa), TRBDF2(nlsolve = nsa))
        JAC_CALLS[] = 0
        sol = solve(prob, alg; reltol = 1.0e-8, abstol = 1.0e-10)
        @test SciMLBase.successful_retcode(sol)
        @test norm(sol.u[end] .- refsol.u[end]) / norm(refsol.u[end]) < 1.0e-4
        # W is reassembled on dt changes throughout (stats.nw ~ O(100) here), but a
        # fresh Jacobian evaluation is only needed on first use and after divergences.
        # Before the J/W split every W update called the user jac, so this tracked
        # stats.nw one-to-one.
        @test JAC_CALLS[] < sol.stats.nw / 3
        @test JAC_CALLS[] >= 1
    end
end

@testset "stale-W Newton directions are rescaled" begin
    # A stage that reuses `W` from `h₀ = W_γdt` at `h = γdt = r h₀` contracts each mode of a
    # linear problem by `1 - s (λ - 1/h)/(λ - 1/h₀)`. The balanced `s = 2r/(1 + r)` gives
    # `(r - 1)/(r + 1) ⋅ (1 + h₀λ)/(1 - h₀λ)`; no correction (`s = 1`) or the reciprocal
    # orientation `s = 2/(1 + r)` give different contractions in at least one regime.
    # The iterates are read off the points `f` is evaluated at, on Prothero-Robinson
    # `u' = λ(u - cos t) - sin t`, whose implicit Euler stage root is known exactly.
    pts = Float64[]
    λ = Ref(0.0)
    function pr!(du, u, p, t)
        u[1] isa Float64 && push!(pts, u[1])
        du[1] = λ[] * (u[1] - cos(t)) - sin(t)
        return nothing
    end
    for nl in (nsa, NLNewton()), λval in (-1.0e6, -1.0e-2), r in (0.9, 1.1)
        λ[] = λval
        integ = init(
            ODEProblem(pr!, [1.0], (0.0, 10.0)), ImplicitEuler(nlsolve = nl);
            dt = 1.0e-2, reltol = 1.0e-2, abstol = 1.0e-2
        )
        for _ in 1:3
            step!(integ)
        end
        h₀ = integ.cache.nlsolver.cache.W_γdt
        set_proposed_dt!(integ, r * h₀)
        t₀, uprev = integ.t, integ.u[1]
        empty!(pts)
        step!(integ)
        h = integ.t - t₀
        iterates = pts[[true; diff(pts) .!= 0]]
        ustar = (uprev - h * (sin(integ.t) + λval * cos(integ.t))) / (1 - h * λval)
        e = iterates .- ustar
        rh = h / h₀
        expected = (rh - 1) / (rh + 1) * (1 + h₀ * λval) / (1 - h₀ * λval)
        @testset "$(nameof(typeof(nl))) λ=$λval r=$r" begin
            @test integ.cache.nlsolver.cache.W_γdt == h₀
            @test rh ≈ r
            @test length(e) >= 2
            @test e[2] / e[1] ≈ expected rtol = 1.0e-4
        end
    end
end

@testset "globalized inner solver does not converge on a rejected step" begin
    # A TrustRegion trial step that is rejected leaves the iterate exactly unmoved
    # without terminating the inner cache; the outer displacement test alone reads
    # that as convergence and accepts the stage with no correction applied (#3817).
    # Driven at a dt far too coarse for the Robertson transient, so the inner solves
    # genuinely cannot converge and the failure must be reported rather than
    # silently absorbed.
    nsa_tr = NonlinearSolveAlg(TrustRegion(; autodiff = AutoForwardDiff()))
    nsa_nr = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff()))
    hard = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1.0e3), [0.04, 3.0e7, 1.0e4])

    sol_tr = solve(hard, FBDF(nlsolve = nsa_tr); dt = 1.0, adaptive = false)
    sol_nr = solve(hard, FBDF(nlsolve = nsa_nr); dt = 1.0, adaptive = false)

    # The state must not be reported as a successful solve while frozen at u0.
    @test !(SciMLBase.successful_retcode(sol_tr) && sol_tr.u[end] == hard.u0)
    # TrustRegion must reach the same verdict as the non-globalized inner solver.
    @test SciMLBase.successful_retcode(sol_tr) == SciMLBase.successful_retcode(sol_nr)
end

@testset "TrustRegion matches NewtonRaphson when the solves do converge" begin
    nsa_tr = NonlinearSolveAlg(TrustRegion(; autodiff = AutoForwardDiff()))
    sol = solve(prob, FBDF(nlsolve = nsa_tr); reltol = 1.0e-8, abstol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)
    @test norm(sol.u[end] .- refsol.u[end]) / norm(refsol.u[end]) < 1.0e-4

    # At default tolerances the rejected-step false convergence used to cost three
    # orders of magnitude (relerr 1.7e-2, where NLNewton gives 2.1e-5); with rejected
    # trial steps excluded from the convergence test it lands back at tolerance level
    # (measured 1.8e-5).
    sol_def = solve(prob, FBDF(nlsolve = nsa_tr))
    @test SciMLBase.successful_retcode(sol_def)
    @test norm(sol_def.u[end] .- refsol.u[end]) / norm(refsol.u[end]) < 1.0e-3
end

@testset "non-adaptive solves do not freeze the Jacobian" begin
    # `do_newJW` returns `(true, true)` whenever `adaptive = false`: with a prescribed
    # step sequence there is no error test to catch the step a stale Jacobian degrades
    # and no way to pull it back. `NonlinearSolveAlg` has to make the same call —
    # freezing `J` for the whole solve here would silently diverge from `NLNewton`.
    fixedprob = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1.0), [0.04, 3.0e7, 1.0e4])
    for nl in (NLNewton(), nsa)
        JAC_CALLS[] = 0
        sol = solve(fixedprob, TRBDF2(nlsolve = nl); dt = 1.0e-4, adaptive = false)
        @test SciMLBase.successful_retcode(sol)
        @test JAC_CALLS[] >= sol.stats.naccept
    end
end
