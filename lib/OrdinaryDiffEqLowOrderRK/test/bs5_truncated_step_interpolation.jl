using OrdinaryDiffEqLowOrderRK, DiffEqBase, Test

# The extra interpolation stages of the non-lazy BS5 interpolant live inside
# kshortsize, so a step shortened by a ContinuousCallback must explicitly force
# them to be recomputed rather than reuse the ones built for the original dt.

exp_prob = ODEProblem((du, u, p, t) -> (du[1] = u[1]; nothing), [1.0], (0.0, 5.0))
root_cb = ContinuousCallback((u, t, integ) -> u[1] - 50.0, integ -> nothing)

function worst_interior_error(sol, lo, hi)
    return maximum(
        abs(sol(t)[1] - exp(t)) / exp(t)
            for t in range(lo, hi, length = 21)[2:(end - 1)]
    )
end

for alg in (BS5(lazy = false), BS5())
    sol = solve(
        exp_prob, alg; abstol = 1.0e-10, reltol = 1.0e-10, callback = root_cb
    )
    i = findfirst(≈(log(50.0)), sol.t)
    @test worst_interior_error(sol, sol.t[i - 1], sol.t[i]) < 1.0e-8
end
