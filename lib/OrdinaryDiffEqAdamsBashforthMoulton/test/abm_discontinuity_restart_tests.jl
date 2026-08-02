using OrdinaryDiffEqAdamsBashforthMoulton, DiffEqBase, Test

# A callback that jumps `u` invalidates the multistep history: the stored `k`
# values describe the pre-jump trajectory, so continuing the Adams formula across
# the jump extrapolates through the discontinuity and excites the parasitic roots.
# The caches guard against this by resetting their step counter whenever
# `integrator.derivative_discontinuity` is set, which restarts the low-order
# startup procedure from the post-jump state.

f_iip(du, u, p, t) = (du[1] = -u[1])
f_oop(u, p, t) = -u

const u0_iip = [10.0]
const u0_oop = 10.0
const tjump = 4.0

reset_iip!(integrator) = (integrator.u[1] = u0_iip[1])
reset_oop!(integrator) = (integrator.u = u0_oop)

# u' = -u restarted from u0 at tjump repeats the initial value problem verbatim, so
# a correctly restarted multistep method must reproduce its own opening segment
# step for step. Any surviving pre-jump history shows up as a mismatch here.
@testset "fixed-step restart: $(nameof(typeof(alg)))" for alg in (
        AB3(), AB4(), AB5(), ABM32(), ABM43(), ABM54(),
    )
    @testset "$name" for (name, f, u0, affect!) in (
            ("in-place", f_iip, u0_iip, reset_iip!),
            ("out-of-place", f_oop, u0_oop, reset_oop!),
        )
        dt = 0.5
        cb = DiscreteCallback((u, t, integrator) -> t == tjump, affect!)
        reference = solve(
            ODEProblem(f, u0, (0.0, tjump)), alg; adaptive = false, dt = dt
        )
        sol = solve(
            ODEProblem(f, u0, (0.0, 2tjump)), alg;
            callback = cb, tstops = [tjump], adaptive = false, dt = dt
        )

        # `save_positions` stores the jump twice; the later save holds the reset state.
        jump = findlast(==(tjump), sol.t)
        @test sol.u[jump] == u0
        for k in 1:(length(reference.t) - 1)
            @test sol.t[jump + k] ≈ tjump + reference.t[1 + k]
            @test sol.u[jump + k] == reference.u[1 + k]
        end
    end
end

# The adaptive variable-coefficient methods carry a variable-order divided-difference
# history rather than a step counter, so check the user-visible consequence instead:
# after the jump the solution must track the exact decay to the requested tolerance.
@testset "adaptive restart: $(nameof(typeof(alg)))" for alg in (
        VCAB3(), VCAB4(), VCAB5(), VCABM3(), VCABM4(), VCABM5(), VCABM(),
    )
    cb = DiscreteCallback((u, t, integrator) -> t == tjump, reset_iip!)
    sol = solve(
        ODEProblem(f_iip, u0_iip, (0.0, 2tjump)), alg;
        callback = cb, tstops = [tjump], abstol = 1.0e-8, reltol = 1.0e-8
    )
    exact(t) = u0_iip[1] * exp(-(t < tjump ? t : t - tjump))
    @test all(
        isapprox(sol.u[i][1], exact(sol.t[i]); atol = 1.0e-5)
            for i in eachindex(sol.t) if sol.t[i] != tjump
    )
end
