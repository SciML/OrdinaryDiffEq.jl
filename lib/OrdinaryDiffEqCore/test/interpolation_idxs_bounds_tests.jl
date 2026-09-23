using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, Test

# In-place dense/continuous interpolation dispatches to per-algorithm
# `_ode_interpolant!` kernels that index `out`/`y₀`/`k` under `@inbounds`. The
# Core entry points must reject mismatched `out`/`idxs` and out-of-range `idxs`
# instead of silently reading or writing out of bounds.
function growth!(du, u, p, t)
    du[1] = u[1]
    du[2] = 2 * u[2]
    return du[3] = -u[3]
end

@testset "in-place interpolation idxs/out validation" begin
    sol = solve(ODEProblem(growth!, [1.0, 2.0, 3.0], (0.0, 1.0)), Tsit5())

    out = zeros(2)
    @test_throws DimensionMismatch sol(out, 0.5; idxs = 1:3)
    @test_throws BoundsError sol(out, 0.5; idxs = [1, 10^6])

    # vector-of-times in-place path
    outs = [zeros(2) for _ in 1:2]
    @test_throws DimensionMismatch sol(outs, [0.4, 0.5]; idxs = 1:3)
    outs_ok = [zeros(2) for _ in 1:2]
    @test_throws BoundsError sol(outs_ok, [0.4, 0.5]; idxs = [1, 10^6])

    # valid `idxs` agrees with the allocating path, which already bounds-checks
    out_ok = zeros(2)
    sol(out_ok, 0.5; idxs = 1:2)
    @test out_ok == sol(0.5; idxs = 1:2)
    @test_throws BoundsError sol(0.5; idxs = [1, 10^6])

    # scalar `idxs` keeps working
    out_scalar = zeros(1)
    sol(out_scalar, 0.5; idxs = 2)
    @test out_scalar[] == sol(0.5; idxs = 2)

    # `idxs = nothing` still writes the full state
    out_full = zeros(3)
    sol(out_full, 0.5)
    @test out_full == sol(0.5)
end
