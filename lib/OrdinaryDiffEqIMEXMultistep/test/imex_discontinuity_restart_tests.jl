using OrdinaryDiffEqIMEXMultistep, DiffEqBase, Test

# CNAB2 and CNLF2 are two-step: they combine the current state with `k2`/`uprev2`
# carried over from the previous step. A callback that jumps `u` invalidates that
# history, so the step after the jump must fall back to the one-step startup
# formula. Before the `is_startup_step` guard, the startup branch was selected by
# `integrator.iter == 1` alone, so the pre-jump history was reused: CNLF2's
# leapfrog step ran from the stale `uprev2` and immediately produced a negative
# value on a problem whose exact solution is a positive decay.

# u' = -u, split evenly between the implicitly and explicitly treated parts.
f1_iip(du, u, p, t) = (du[1] = -0.5 * u[1])
f2_iip(du, u, p, t) = (du[1] = -0.5 * u[1])
f1_oop(u, p, t) = -0.5 * u
f2_oop(u, p, t) = -0.5 * u

const u0_iip = [10.0]
const u0_oop = 10.0
const tjump = 4.0

reset_iip!(integrator) = (integrator.u[1] = u0_iip[1])
reset_oop!(integrator) = (integrator.u = u0_oop)

# Resetting to u0 at tjump repeats the original initial value problem, so a
# correctly restarted two-step method reproduces its own opening segment exactly.
@testset "$(nameof(typeof(alg)))" for alg in (CNAB2(), CNLF2())
    @testset "$name" for (name, f1, f2, u0, affect!) in (
            ("in-place", f1_iip, f2_iip, u0_iip, reset_iip!),
            ("out-of-place", f1_oop, f2_oop, u0_oop, reset_oop!),
        )
        dt = 0.5
        cb = DiscreteCallback((u, t, integrator) -> t == tjump, affect!)
        reference = solve(
            SplitODEProblem(f1, f2, u0, (0.0, tjump)), alg; adaptive = false, dt
        )
        sol = solve(
            SplitODEProblem(f1, f2, u0, (0.0, 2tjump)), alg;
            callback = cb, tstops = [tjump], adaptive = false, dt
        )

        # `save_positions` stores the jump twice; the later save holds the reset state.
        jump = findlast(==(tjump), sol.t)
        @test sol.u[jump] == u0
        for k in 1:(length(reference.t) - 1)
            @test sol.t[jump + k] ≈ tjump + reference.t[1 + k]
            @test sol.u[jump + k] == reference.u[1 + k]
        end

        # The exact solution is a positive decay on both sides of the jump.
        @test all(all(>(0), u) for u in sol.u)
    end
end
