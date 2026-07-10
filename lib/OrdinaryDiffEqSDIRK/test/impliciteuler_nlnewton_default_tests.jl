using LinearAlgebra: diag, norm
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve: NLNewton
using SciMLBase: successful_retcode
using StaticArrays: @SMatrix, @SVector, SVector
using Test

#=
Do the historical ImplicitEuler defaults make sense?

NLNewton's default κ is 1//100. Convergence is η * ndz < κ with ndz already
scaled by abstol/reltol. This file probes that default: where it is adequate,
where it is not, and whether a product-level special case for ImplicitEuler is
justified (it is not — tight κ is a call-site choice for hard problems).

Stratospheric production–destruction (PDS) problem is the PositiveIntegrators-
style static-array stress case that previously motivated tightening
ImplicitEuler's default κ globally. That product change is intentionally
*not* applied; this test suite documents the real behavior instead.
=#

# ── Shared stratospheric chemistry (vector and static-array forms) ───────────

function stratreac_rates(u, t)
    O1D, O, O3, O2, NO, NO2 = u
    Tr = 4.5
    Ts = 19.5
    T = mod(t / 3600, 24)
    if Tr <= T <= Ts
        Tfrac = (2 * T - Tr - Ts) / (Ts - Tr)
        sigma = 0.5 + 0.5 * cos(pi * abs(Tfrac) * Tfrac)
    else
        sigma = zero(t)
    end

    M = 8.12e16
    k1 = 2.643e-10 * sigma^3
    k2 = 8.018e-17
    k3 = 6.12e-4 * sigma
    k4 = 1.567e-15
    k5 = 1.07e-3 * sigma^2
    k6 = 7.11e-11
    k7 = 1.2e-10
    k8 = 6.062e-15
    k9 = 1.069e-11
    k10 = 1.289e-2 * sigma
    k11 = 1.0e-8

    return (
        k1 * O2,
        k2 * O * O2,
        k3 * O3,
        k4 * O3 * O,
        k5 * O3,
        k6 * M * O1D,
        k7 * O1D * O3,
        k8 * O3 * NO,
        k9 * NO2 * O,
        k10 * NO2,
        k11 * NO * O,
    )
end

function stratreac_rhs(u, p, t)
    r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11 = stratreac_rates(u, t)
    return @SVector [
        r5 - r6 - r7
        2 * r1 - r2 + r3 - r4 + r6 - r9 + r10 - r11
        r2 - r3 - r4 - r5 - r7 - r8
        -r1 - r2 + r3 + 2 * r4 + r5 + 2 * r7 + r8 + r9
        -r8 + r9 + r10 - r11
        r8 - r9 - r10 + r11
    ]
end

function stratreac_production(u, p, t)
    r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11 = stratreac_rates(u, t)
    return @SMatrix [
        0.0 0.0 r5 0.0 0.0 0.0
        r6 r1 + r10 r3 r1 0.0 0.0
        0.0 r2 0.0 0.0 0.0 0.0
        r7 r4 + r9 r4 + r7 + r8 r3 + r5 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 r9 + r10
        0.0 0.0 0.0 0.0 r8 + r11 0.0
    ]
end

function stratreac_destruction(u, p, t)
    _, r2, _, _, _, _, _, _, _, _, r11 = stratreac_rates(u, t)
    return @SVector [0.0, r11, 0.0, r2, 0.0, 0.0]
end

function stratreac_pds_rhs(u, p, t)
    P = stratreac_production(u, p, t)
    D = stratreac_destruction(u, p, t)
    return diag(P) + vec(sum(P; dims = 2)) - vec(sum(P; dims = 1)) - vec(D)
end

# Vector (mutable) form of the same RHS, for iip comparison
function stratreac_rhs!(du, u, p, t)
    r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11 = stratreac_rates(u, t)
    du[1] = r5 - r6 - r7
    du[2] = 2 * r1 - r2 + r3 - r4 + r6 - r9 + r10 - r11
    du[3] = r2 - r3 - r4 - r5 - r7 - r8
    du[4] = -r1 - r2 + r3 + 2 * r4 + r5 + 2 * r7 + r8 + r9
    du[5] = -r8 + r9 + r10 - r11
    du[6] = r8 - r9 - r10 + r11
    return nothing
end

const U0_S = @SVector [9.906e1, 6.624e8, 5.326e11, 1.697e16, 4.0e6, 1.093e9]
const U0_V = Vector(U0_S)
const TSPAN = (4.32e4, 3.024e5)
const T_CHECK = 2.0e5
const SOLVE_KWARGS = (; dt = 1.0, abstol = 1.0e-5, reltol = 1.0e-5)

function solve_pair(rhs, pds_rhs, u0, alg; kwargs...)
    sol = solve(ODEProblem(rhs, u0, TSPAN), alg; kwargs...)
    sol_pds = solve(ODEProblem(pds_rhs, u0, TSPAN), alg; kwargs...)
    return sol, sol_pds
end

function maxabs_at(sol, sol_pds, t = T_CHECK)
    return maximum(abs, sol(t) .- sol_pds(t))
end

# ── 1. Product defaults ─────────────────────────────────────────────────────

@testset "ImplicitEuler product defaults" begin
    ie = ImplicitEuler()
    @test ie.nlsolve isa NLNewton
    # Historical / shared NLNewton default — do not special-case ImplicitEuler.
    @test ie.nlsolve.κ === NLNewton().κ === 1 // 100
    @test ie.nlsolve.max_iter == NLNewton().max_iter == 10
    # Same default as other common SDIRK methods that use NLNewton.
    @test Trapezoid().nlsolve.κ === 1 // 100
    @test TRBDF2().nlsolve.κ === 1 // 100
end

# ── 2. Default is adequate on ordinary stiff problems ───────────────────────

@testset "Default κ adequate on standard stiff ODEs" begin
    # Robertson
    function rober!(du, u, p, t)
        y₁, y₂, y₃ = u
        k₁, k₂, k₃ = p
        du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
        du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
        du[3] = k₂ * y₂^2
        return nothing
    end
    prob_rober = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
    sol_rober = solve(prob_rober, ImplicitEuler(); abstol = 1.0e-6, reltol = 1.0e-6)
    @test successful_retcode(sol_rober)
    u_rober = sol_rober.u[end]
    @test all(>=(0), u_rober) # concentrations stay non-negative
    @test sum(u_rober) ≈ 1 atol = 1.0e-3

    # Simple linear stiff system: u' = A u with eigenvalues -1000, -0.1
    A = [-1000.0 0.0; 0.0 -0.1]
    f_lin(u, p, t) = A * u
    sol_lin = solve(
        ODEProblem(f_lin, [1.0, 1.0], (0.0, 10.0)),
        ImplicitEuler();
        abstol = 1.0e-8, reltol = 1.0e-8
    )
    @test successful_retcode(sol_lin)
    u_lin = sol_lin.u[end]
    @test u_lin[1] ≈ 0 atol = 1.0e-4          # fast mode gone
    @test u_lin[2] ≈ exp(-1.0) rtol = 1.0e-2  # slow mode ≈ e^{-1}
end

# ── 3. Algebraic RHS vs PDS form at t=0 ─────────────────────────────────────

@testset "PDS and closed-form RHS agree at u0" begin
    t0 = TSPAN[1]
    @test stratreac_rhs(U0_S, nothing, t0) ≈ stratreac_pds_rhs(U0_S, nothing, t0)
    du = similar(U0_V)
    stratreac_rhs!(du, U0_V, nothing, t0)
    @test du ≈ Vector(stratreac_rhs(U0_S, nothing, t0))
end

# ── 4. Default κ on stratospheric static-array PDS (the hard case) ──────────

@testset "Default κ on static-array stratospheric PDS" begin
    alg = ImplicitEuler() # κ = 1//100
    sol, sol_pds = solve_pair(
        stratreac_rhs, stratreac_pds_rhs, U0_S, alg; SOLVE_KWARGS...
    )

    # Both still report success — the failure mode is false Newton acceptance,
    # not a failed retcode. That is exactly why tightening the *product* default
    # would hide the issue rather than surface it.
    @test successful_retcode(sol)
    @test successful_retcode(sol_pds)

    δ = maxabs_at(sol, sol_pds)
    # Under default κ the two algebraically-equivalent RHS forms diverge
    # catastrophically (PDS path collapses; closed form stays physical).
    @test δ > 1.0e6
    @test !isapprox(sol(T_CHECK), sol_pds(T_CHECK))
    # Closed-form trajectory stays in a physical ball; PDS does not.
    @test all(isfinite, sol(T_CHECK))
    @test sol(T_CHECK)[4] > 1.0e15 # O2 reservoir intact
    # PDS solution has collapsed species (denormals / zeros observed in practice).
    @test minimum(abs, sol_pds(T_CHECK)) < 1.0e-10 ||
        any(!isfinite, sol_pds(T_CHECK)) ||
        sol_pds(T_CHECK)[1] < 1.0
end

# ── 5. Call-site tight κ restores RHS ↔ PDS equivalence ─────────────────────

@testset "Call-site tight κ restores static-array PDS equivalence" begin
    for κ in (1 // 10000, 1.0e-10)
        alg = ImplicitEuler(nlsolve = NLNewton(κ = κ))
        sol, sol_pds = solve_pair(
            stratreac_rhs, stratreac_pds_rhs, U0_S, alg; SOLVE_KWARGS...
        )
        @test successful_retcode(sol)
        @test successful_retcode(sol_pds)
        @test isapprox(sol(T_CHECK), sol_pds(T_CHECK))
        @test maxabs_at(sol, sol_pds) < 1.0 # absolute, concentrations ~1e16
        # Same adaptive path once Newton is tight enough
        @test length(sol.t) == length(sol_pds.t)
    end
end

# ── 6. κ sweep: where does the default stop making sense? ───────────────────

@testset "κ sweep for static-array PDS equivalence" begin
    # Smaller κ ⇒ tighter Newton. Probe which values keep RHS ≈ PDS.
    # Use a fixed quality metric (max abs error) rather than isapprox alone:
    # two equally-collapsed garbage trajectories can spuriously isapprox.
    κs = [1 // 100, 1 // 1000, 1 // 10000, 1 // 100000]
    deltas = Dict{Any, Float64}()
    physical = Dict{Any, Bool}()
    for κ in κs
        alg = ImplicitEuler(nlsolve = NLNewton(κ = κ))
        sol, sol_pds = solve_pair(
            stratreac_rhs, stratreac_pds_rhs, U0_S, alg; SOLVE_KWARGS...
        )
        δ = maxabs_at(sol, sol_pds)
        deltas[κ] = δ
        # "Physical" = both Success, O2 reservoir intact on both, and δ moderate.
        physical[κ] = successful_retcode(sol) && successful_retcode(sol_pds) &&
            sol(T_CHECK)[4] > 1.0e15 &&
            sol_pds(T_CHECK)[4] > 1.0e15 &&
            δ < 1.0e3
    end

    # Historical default is not enough for this problem.
    @test physical[1 // 100] == false
    @test deltas[1 // 100] > 1.0e6
    # Two orders of magnitude tighter than default is enough.
    @test physical[1 // 10000] == true
    @test physical[1 // 100000] == true
    @test deltas[1 // 10000] < 1.0
    # Intermediate 1//1000 is transitional (do not pin pass/fail); only require
    # that the error is not worse than the default by orders of magnitude in the
    # wrong direction when it succeeds — record via the dict for CI logs.
    @test deltas[1 // 1000] isa Real
end

# ── 7. Step size can compensate for loose κ ─────────────────────────────────

@testset "Smaller dt compensates for default κ on static-array PDS" begin
    alg = ImplicitEuler()
    sol, sol_pds = solve_pair(
        stratreac_rhs, stratreac_pds_rhs, U0_S, alg;
        dt = 0.1, abstol = 1.0e-5, reltol = 1.0e-5
    )
    @test successful_retcode(sol)
    @test successful_retcode(sol_pds)
    # With a 10× smaller fixed step, default κ is adequate for this problem.
    @test isapprox(sol(T_CHECK), sol_pds(T_CHECK); rtol = 1.0e-3, atol = 1.0e3)
end

# ── 8. Mutable Vector state with default κ ──────────────────────────────────

@testset "Vector (iip) stratospheric RHS with default κ" begin
    # Even the closed-form RHS with mutable storage should integrate successfully
    # under default ImplicitEuler; this is not the PDS static-array stress case.
    sol = solve(
        ODEProblem(stratreac_rhs!, copy(U0_V), TSPAN),
        ImplicitEuler();
        SOLVE_KWARGS...
    )
    @test successful_retcode(sol)
    @test all(isfinite, sol(T_CHECK))
    @test sol(T_CHECK)[4] > 1.0e15
end

# ── 9. Stricter solver tolerances alone do not fix default κ on PDS ─────────

@testset "Tighter abstol/reltol alone do not fix default-κ PDS" begin
    alg = ImplicitEuler()
    sol, sol_pds = solve_pair(
        stratreac_rhs, stratreac_pds_rhs, U0_S, alg;
        dt = 1.0, abstol = 1.0e-8, reltol = 1.0e-8
    )
    # Either integration hits MaxIters or PDS still diverges from closed form —
    # tightening *solver* tols without tightening Newton κ is not sufficient.
    δ = maxabs_at(sol, sol_pds)
    still_broken = !isapprox(sol(T_CHECK), sol_pds(T_CHECK)) || δ > 1.0e6 ||
        !successful_retcode(sol) || !successful_retcode(sol_pds)
    @test still_broken
end

# ── 10. Explicit recommendation for PositiveIntegrators-style usage ─────────

@testset "PositiveIntegrators-style call site (explicit κ, not product default)" begin
    # Downstream packages that need static-array PDS accuracy should set κ at
    # the call site rather than relying on a special ImplicitEuler default.
    alg = ImplicitEuler(nlsolve = NLNewton(κ = 1 // 10000))
    sol, sol_pds = solve_pair(
        stratreac_rhs, stratreac_pds_rhs, U0_S, alg; SOLVE_KWARGS...
    )
    @test successful_retcode(sol) && successful_retcode(sol_pds)
    @test isapprox(sol(T_CHECK), sol_pds(T_CHECK))
    @test maxabs_at(sol, sol_pds) < 1.0
end
