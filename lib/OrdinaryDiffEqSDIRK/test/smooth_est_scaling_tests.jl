using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqCore: get_EEst
using SciMLBase
using Test

# Issue #2902: smooth_est must solve W ERR = W_γdt⁻¹ err (Hairer–Wanner / Shampine),
# using the γdt stored with the (possibly reused) W factorization.
#
# Non-stiff (J ≡ 0): W = −W_γdt⁻¹ I ⇒ ‖ERR‖ = ‖err‖ (ratio ≈ 1).
# Stiff linear (J = −λ): M = I + W_γdt·λ ⇒ ‖ERR‖/‖err‖ = 1/(1 + λ·W_γdt).

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
            # Pre-fix buggy ratio was (γ·dt)/(1 + λ·γ·dt) = (γ·dt)·expected
            buggy = Float64(γ) * dt * expected
            @test abs(ratio - buggy) > 1.0e-3 * max(abs(expected), abs(buggy))
        end
    end
end

# Compare the smoothed estimate to the true local error u_num - u(t₀+dt) on u' = λu.
# With abstol=1, reltol=0 the residual weight is 1, so EEst equals ‖est‖.
# Non-stiff (z = λ·dt ∈ {-0.1, -1}): TRBDF2 ratio ≈ 1.0 / 1.07 — stay in a derived band.
# Stiff (λ = -1e6): smoothed must be filtered well below the raw estimate
# (exact exp(λ·dt) underflows in Float64, so no true-error ratio there).
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
    # abstol=1, reltol=0 ⇒ EEst = |est| for a scalar problem
    tol_kwargs = (; abstol = 1.0, reltol = 0.0)

    @testset "non-stiff band (TRBDF2)" begin
        # Reviewer-measured TRBDF2 ratios ≈ 1.0 at z=-0.1 and ≈ 1.07 at z=-1.
        for (z, lo, hi) in ((-0.1, 0.85, 1.15), (-1.0, 0.9, 1.25))
            λ = z # with dt=1 ⇒ z = λ·dt
            dt = 1.0
            prob = ODEProblem(ODEFunction(f_lin!; jac = jac_lin!), u0, (0.0, dt), [λ])
            es, integ = first_step_EEst(
                prob, TRBDF2(smooth_est = true); dt = dt, tol_kwargs...
            )
            true_err = abs(integ.u[1] - u0[1] * exp(λ * dt))
            @test true_err > 0
            ratio = es / true_err
            @test lo <= ratio <= hi
        end
    end

    @testset "stiff filtering (λ = -1e6)" begin
        λ = -1.0e6
        dt = 1.0e-2
        prob = ODEProblem(ODEFunction(f_lin!; jac = jac_lin!), u0, (0.0, 1.0), [λ])
        for ALG in (TRBDF2, KenCarp4, SDIRK2)
            es, _ = first_step_EEst(
                prob, ALG(smooth_est = true); dt = dt, tol_kwargs...
            )
            er, _ = first_step_EEst(
                prob, ALG(smooth_est = false); dt = dt, tol_kwargs...
            )
            @test isfinite(es) && isfinite(er) && es > 0 && er > 0
            # Smoothing must filter the explicit-like raw estimate
            @test es < er / 10
        end
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
