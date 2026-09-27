using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqCore: get_EEst
using SciMLBase
using Test

# Issue #2902: smooth_est must solve W ERR = (hγ)⁻¹ err, not W ERR = err.
#
# Non-stiff (J ≡ 0): W = −(hγ)⁻¹ I, so the correctly scaled solve gives ‖ERR‖ = ‖err‖
# (ratio ≈ 1). The buggy W\err path instead yields ratio = γ·dt.
#
# Stiff linear (J = −λ): M = I + hγλ, so ‖ERR‖/‖err‖ = 1/(1 + λ·γ·dt).

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

function first_step_EEst(prob, alg; dt)
    integ = init(
        prob, alg; dt = dt, adaptive = true, dtmax = dt, dtmin = dt, force_dtmin = true
    )
    step!(integ)
    return get_EEst(integ)
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
                es = first_step_EEst(prob, ALG(smooth_est = true); dt = dt)
                er = first_step_EEst(prob, ALG(smooth_est = false); dt = dt)
                @test isfinite(es) && isfinite(er) && er > 0
                ratio = es / er
                # Correct Hairer–Wanner scaling: smoothed ≡ raw when J=0
                @test ratio ≈ 1 atol = 1.0e-8 rtol = 1.0e-6
                @info "$name $label dt=$dt" ratio = ratio pre_fix_γdt = Float64(γ) * dt
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
            es = first_step_EEst(prob, ALG(smooth_est = true); dt = dt)
            er = first_step_EEst(prob, ALG(smooth_est = false); dt = dt)
            @test isfinite(es) && isfinite(er) && er > 0
            ratio = es / er
            @test ratio ≈ expected rtol = 1.0e-6 atol = 1.0e-10
            # Pre-fix buggy ratio was (γ·dt)/(1 + λ·γ·dt) = (γ·dt)·expected
            buggy = Float64(γ) * dt * expected
            @test abs(ratio - buggy) > 1.0e-3 * max(abs(expected), abs(buggy))
            @info "$name $label stiff" ratio = ratio expected = expected
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
        @info "$ALG stiff steps" smooth = s_smooth.stats.naccept raw = s_raw.stats.naccept
    end
end
