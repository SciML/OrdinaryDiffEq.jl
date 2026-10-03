using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqCore: get_EEst
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit
using LinearAlgebra
using SciMLBase
using Test

# Issue #2902: smooth_est must solve W ERR = W_γdt⁻¹ (M err) (Hairer–Wanner / Shampine),
# using the γdt stored with the (possibly reused) W factorization.
#
# Non-stiff (J ≡ 0, M = I): W = −W_γdt⁻¹ I ⇒ ‖ERR‖ = ‖err‖ (ratio ≈ 1).
# Stiff linear (J = −λ, M = I): ‖ERR‖/‖err‖ = 1/(1 + λ·W_γdt).

function f_cos!(du, u, p, t)
    du[1] = cos(t)
    return nothing
end
function jac0!(J, u, p, t)
    J[1, 1] = 0.0
    return nothing
end
f_cos(u, p, t) = [cos(t)]
jac0(u, p, t) = zeros(1, 1)

function first_step_EEst(prob, alg; dt, kwargs...)
    integ = init(
        prob, alg; dt = dt, adaptive = true, dtmax = dt, dtmin = dt, force_dtmin = true,
        kwargs...
    )
    step!(integ)
    return get_EEst(integ), integ
end

@testset "smooth_est scaling with J=0 (#2902)" begin
    # γ = A_{s,s} from each method's tableau; pre-fix ratio was exactly γ·dt.
    algs = (
        ("KenCarp4", KenCarp4, 1 // 4),
        ("KenCarp47", KenCarp47, 1235 // 10000),
        ("TRBDF2", TRBDF2, 1 - 1 / sqrt(2)),
        ("SDIRK2", SDIRK2, 1),
        ("Kvaerno5", Kvaerno5, 0.26),
    )
    for (name, ALG, γ) in algs
        for (label, prob) in (
                ("iip", ODEProblem(ODEFunction(f_cos!; jac = jac0!), [0.0], (0.0, 1.0))),
                ("oop", ODEProblem(ODEFunction(f_cos; jac = jac0), [0.0], (0.0, 1.0))),
            )
            for dt in (0.1, 0.01)
                es, _ = first_step_EEst(prob, ALG(smooth_est = true); dt = dt)
                er, _ = first_step_EEst(prob, ALG(smooth_est = false); dt = dt)
                @test isfinite(es) && isfinite(er) && er > 0
                ratio = es / er
                # Correct Hairer–Wanner scaling: smoothed ≡ raw when J=0
                @test ratio ≈ 1 atol = 1.0e-8 rtol = 1.0e-6
            end
        end
    end
end

@testset "smooth_est stiff linear scaling (#2902)" begin
    λ = 1.0e6
    function f_stiff!(du, u, p, t)
        du[1] = -p[1] * u[1]
        return nothing
    end
    function jac_stiff!(J, u, p, t)
        J[1, 1] = -p[1]
        return nothing
    end
    f_stiff(u, p, t) = [-p[1] * u[1]]
    jac_stiff(u, p, t) = fill(-p[1], 1, 1)

    algs = (
        ("KenCarp4", KenCarp4, 1 // 4),
        ("TRBDF2", TRBDF2, 1 - 1 / sqrt(2)),
        ("SDIRK2", SDIRK2, 1),
    )
    dt = 0.01
    for (name, ALG, γ) in algs
        expected = 1 / (1 + λ * Float64(γ) * dt)
        for (label, prob) in (
                (
                    "iip",
                    ODEProblem(
                        ODEFunction(f_stiff!; jac = jac_stiff!), [1.0], (0.0, 1.0), [λ]
                    ),
                ),
                (
                    "oop",
                    ODEProblem(
                        ODEFunction(f_stiff; jac = jac_stiff), [1.0], (0.0, 1.0), [λ]
                    ),
                ),
            )
            es, _ = first_step_EEst(prob, ALG(smooth_est = true); dt = dt)
            er, _ = first_step_EEst(prob, ALG(smooth_est = false); dt = dt)
            @test isfinite(es) && isfinite(er) && er > 0
            ratio = es / er
            @test ratio ≈ expected rtol = 1.0e-6 atol = 1.0e-10
            buggy = Float64(γ) * dt * expected
            @test abs(ratio - buggy) > 1.0e-3 * max(abs(expected), abs(buggy))
        end
    end
end

# Scalar linear u' = λu, M = I. With abstol=1, reltol=0, EEst = |est|.
# Closed form: W = λ − 1/(γ_diag·dt), W ERR = (γ_diag·dt)⁻¹ err
# ⇒ |ERR| = |err| / |1 − γ_diag·z| with z = λ·dt. True local error = |φ(z) − eᶻ|
# (when eᶻ underflows to 0, that is just |φ(z)| = |u_num|).
@testset "smooth_est vs true local error (#2902)" begin
    function f_lin!(du, u, p, t)
        du[1] = p[1] * u[1]
        return nothing
    end
    function jac_lin!(J, u, p, t)
        J[1, 1] = p[1]
        return nothing
    end
    u0 = [1.0]
    tol_kwargs = (; abstol = 1.0, reltol = 0.0)
    # TRBDF2 diagonal entry γ_diag = A_{s,s} = 1 − 1/√2
    γ_diag = 1 - 1 / sqrt(2)

    @testset "closed form (TRBDF2)" begin
        for z in (-0.1, -1.0)
            λ = z
            dt = 1.0
            prob = ODEProblem(ODEFunction(f_lin!; jac = jac_lin!), u0, (0.0, dt), [λ])
            es, integ = first_step_EEst(
                prob, TRBDF2(smooth_est = true); dt = dt, tol_kwargs...
            )
            er, _ = first_step_EEst(
                prob, TRBDF2(smooth_est = false); dt = dt, tol_kwargs...
            )
            φ = integ.u[1] / u0[1]
            e_z = exp(z)
            true_err = abs(φ - e_z)
            @test true_err > 0 && er > 0
            # Smoothed estimate from the raw embedded difference
            smooth_closed = er / abs(1 - γ_diag * z)
            # rtol from the order-2 leading-term cancellation: |φ−eᶻ| ~ O(z³) so a few
            # ulps on the estimate are enough once the closed form is exact.
            @test es ≈ smooth_closed rtol = 1.0e-9
            # Distinct from the line above: estimate tracks the true local truncation
            # error to within a small O(1) factor set by the order-2 leading constant
            # (not a restatement of es ≈ smooth_closed).
            @test 0.25 < es / true_err < 4.0
        end
    end

    @testset "stiff filtering (λ = -1e6)" begin
        λ = -1.0e6
        dt = 1.0e-2
        z = λ * dt
        prob = ODEProblem(ODEFunction(f_lin!; jac = jac_lin!), u0, (0.0, 1.0), [λ])
        for (ALG, γd) in ((TRBDF2, γ_diag), (KenCarp4, 1 // 4), (SDIRK2, 1))
            es, integ = first_step_EEst(
                prob, ALG(smooth_est = true); dt = dt, tol_kwargs...
            )
            er, _ = first_step_EEst(
                prob, ALG(smooth_est = false); dt = dt, tol_kwargs...
            )
            # eᶻ underflows to 0 ⇒ true local error is |u_num|
            true_err = abs(integ.u[1] - u0[1] * exp(z))
            @test true_err > 0 && es > 0 && er > 0
            @test es < er / 10
            @test es ≈ er / abs(1 - Float64(γd) * z) rtol = 1.0e-9
        end
    end
end

# Singular mass-matrix Robertson DAE: without M·err premultiply, smooth_est blows up
# on the algebraic row and the solver returns Unstable.
@testset "smooth_est singular mass matrix Robertson (#2902)" begin
    function rober!(du, u, p, t)
        y₁, y₂, y₃ = u
        k₁, k₂, k₃ = p
        du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
        du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
        du[3] = y₁ + y₂ + y₃ - 1
        return nothing
    end
    M = Diagonal([1.0, 1.0, 0.0])
    prob = ODEProblem(
        ODEFunction{true}(rober!; mass_matrix = M),
        [1.0, 0.0, 0.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4)
    )
    ref = solve(
        prob, KenCarp4(smooth_est = false); reltol = 1.0e-12, abstol = 1.0e-12,
        initializealg = BrownFullBasicInit()
    )
    @test SciMLBase.successful_retcode(ref)
    for ALG in (SDIRK2, Hairer4, Hairer42, KenCarp4, Kvaerno5, TRBDF2)
        sol = solve(
            prob, ALG(smooth_est = true); reltol = 1.0e-8, abstol = 1.0e-8,
            initializealg = BrownFullBasicInit(), maxiters = 100_000
        )
        @test sol.retcode == ReturnCode.Success
        # atol from the solve abstol so the small y₂ (~1e-7) component is checked;
        # rtol covers the O(1) components (y₁, y₃).
        @test sol.u[end] ≈ ref(sol.t[end]) rtol = 1.0e-3 atol = 1.0e-8
    end
end

# Matrix-shaped u with a non-identity diagonal mass matrix: M·err must go through
# _vec or mul! hits DimensionMismatch on the (length(u)×length(u)) mass matrix.
# Use implicit-first-stage methods (SDIRK2/Hairer4): explicit-first-stage tableaus
# still hit a separate `_mmdiv` broadcast issue with matrix u (unrelated to #2902).
@testset "smooth_est matrix state with mass matrix (#2902)" begin
    function f_mat!(du, u, p, t)
        @. du = -u
        return nothing
    end
    u0 = ones(2, 2)
    M = Diagonal([1.0, 1.0, 1.0, 0.5])
    prob = ODEProblem(ODEFunction(f_mat!; mass_matrix = M), u0, (0.0, 0.1))
    for ALG in (SDIRK2, Hairer4)
        sol = solve(
            prob, ALG(smooth_est = true); reltol = 1.0e-6, abstol = 1.0e-6,
            dt = 1.0e-2
        )
        @test SciMLBase.successful_retcode(sol)
        @test size(sol.u[end]) == size(u0)
    end
end

# Stiff Van der Pol: raw embedded estimate is explicit-like; correctly scaled smoothing
# must change adaptive behaviour (accept fewer steps than the raw estimate).
@testset "smooth_est stiff step counts (#2902)" begin
    function vdp!(du, u, p, t)
        du[1] = u[2]
        du[2] = p[1] * ((1 - u[1]^2) * u[2] - u[1])
        return nothing
    end
    prob = ODEProblem(vdp!, [2.0, 0.0], (0.0, 6.3), [1.0e5])
    for ALG in (TRBDF2, KenCarp4)
        s_smooth = solve(prob, ALG(); reltol = 1.0e-8, abstol = 1.0e-11)
        s_raw = solve(prob, ALG(smooth_est = false); reltol = 1.0e-8, abstol = 1.0e-11)
        @test SciMLBase.successful_retcode(s_smooth)
        @test SciMLBase.successful_retcode(s_raw)
        @test s_smooth.stats.naccept < s_raw.stats.naccept
    end
end
