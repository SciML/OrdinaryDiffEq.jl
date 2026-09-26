using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, RecursiveArrayTools, Test

# In-place interpolation must reject mismatched/out-of-range `idxs`/`out` before `@inbounds` kernels.
function growth!(du, u, p, t)
    du[1] = u[1]
    du[2] = 2 * u[2]
    return du[3] = -u[3]
end

@testset "in-place interpolation idxs/out validation" begin
    prob = ODEProblem(growth!, [1.0, 2.0, 3.0], (0.0, 1.0))
    sol = solve(prob, Tsit5())

    out = zeros(2)
    @test_throws DimensionMismatch sol(out, 0.5; idxs = 1:3)
    @test_throws BoundsError sol(out, 0.5; idxs = [1, 10^6])

    # vector-of-times in-place path, including mixed-length outputs
    outs = [zeros(2) for _ in 1:2]
    @test_throws DimensionMismatch sol(outs, [0.4, 0.5]; idxs = 1:3)
    outs_ok = [zeros(2) for _ in 1:2]
    @test_throws BoundsError sol(outs_ok, [0.4, 0.5]; idxs = [1, 10^6])
    outs_mixed = [zeros(2), zeros(1)]
    @test_throws DimensionMismatch sol(outs_mixed, [0.4, 0.5]; idxs = 1:2)
    outs_voa = VectorOfArray([zeros(2) for _ in 1:3])
    ts3 = [0.3, 0.4, 0.5]
    sol(outs_voa, ts3; idxs = 1:2)
    @test outs_voa.u == sol(ts3; idxs = 1:2).u

    # integrator-path interpolation
    integ = init(prob, Tsit5())
    step!(integ)
    t_mid = integ.t - integ.dt / 2
    @test_throws DimensionMismatch integ(zeros(2), t_mid; idxs = 1:3)
    @test_throws BoundsError integ(zeros(2), t_mid; idxs = [1, 10^6])

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

    # `idxs = nothing` must reject wrong-length `out` (short truncates; long OOB-reads state)
    @test_throws DimensionMismatch sol(zeros(2), 0.5)
    @test_throws DimensionMismatch sol(zeros(4), 0.5)
    @test_throws DimensionMismatch sol([zeros(2), zeros(3)], [0.4, 0.5])
    @test_throws DimensionMismatch sol([zeros(4), zeros(3)], [0.4, 0.5])
    # only a later output has the wrong length
    @test_throws DimensionMismatch sol([zeros(3), zeros(2)], [0.4, 0.5])
    @test_throws DimensionMismatch integ(zeros(2), t_mid)
    @test_throws DimensionMismatch integ(zeros(4), t_mid)

    # extrapolation path (current_extrapolant! / ode_extrapolant!)
    step!(integ)
    t_ex = integ.t + 0.01
    @test_throws DimensionMismatch OrdinaryDiffEqCore.current_extrapolant!(
        zeros(2), t_ex, integ
    )
    @test_throws DimensionMismatch OrdinaryDiffEqCore.current_extrapolant!(
        zeros(4), t_ex, integ
    )
    @test_throws DimensionMismatch OrdinaryDiffEqCore.ode_extrapolant!(
        zeros(2), 1.1, integ, nothing, Val{0}
    )
    @test_throws DimensionMismatch OrdinaryDiffEqCore.ode_extrapolant!(
        zeros(4), 1.1, integ, nothing, Val{0}
    )
end

@testset "in-place interpolation idxs/out validation with Bool masks" begin
    f3!(du, u, p, t) = (du .= u)
    s = solve(ODEProblem(f3!, [1.0, 2.0, 3.0], (0.0, 1.0)), Tsit5(); dense = false)
    for mask in ([true, false, true], BitVector([true, false, true]))
        out_m = zeros(2)
        s(out_m, 0.5; idxs = mask)
        @test out_m == s(0.5; idxs = mask)
        outs_m = [zeros(2)]
        s(outs_m, [0.5]; idxs = mask)
        @test outs_m[1] == only(s([0.5]; idxs = mask).u)
        @test_throws DimensionMismatch s(zeros(3), 0.5; idxs = mask)
        @test_throws DimensionMismatch s([zeros(1)], [0.5]; idxs = mask)
        @test_throws BoundsError s(zeros(3), 0.5; idxs = [mask; true])
    end
end

@testset "in-place interpolation idxs/out validation after resize" begin
    f2!(du, u, p, t) = (du .= u)

    # state grown 2 -> 3 at t = 0.5
    function grow_affect!(i)
        resize!(i, 3)
        i.u[3] = 2.0
        return nothing
    end
    s_grow = solve(
        ODEProblem(f2!, [1.0, 2.0], (0.0, 1.0)), Tsit5();
        callback = DiscreteCallback((u, t, i) -> t == 0.5, grow_affect!),
        tstops = [0.5]
    )
    @test length(s_grow(0.75)) == 3

    out_g = zeros(1)
    s_grow(out_g, 0.75; idxs = [3])
    @test out_g[] == only(s_grow(0.75; idxs = [3]))
    out_g2 = zeros(2)
    s_grow(out_g2, 0.75; idxs = 2:3)
    @test out_g2 == s_grow(0.75; idxs = 2:3)
    @test_throws BoundsError s_grow(zeros(2), 0.75; idxs = [3, 4])
    # before the resize the state still has length 2
    @test_throws BoundsError s_grow(zeros(1), 0.25; idxs = [3])
    # at the resize time, :left sees the pre-resize state, :right the post-resize one
    @test_throws BoundsError s_grow(zeros(1), 0.5; idxs = [3], continuity = :left)
    out_gr = zeros(1)
    s_grow(out_gr, 0.5; idxs = [3], continuity = :right)
    @test out_gr[] == only(s_grow(0.5; idxs = [3], continuity = :right))
    # vector-of-times validates each time's own interval
    @test_throws BoundsError s_grow([zeros(1), zeros(1)], [0.25, 0.75]; idxs = [3])
    outs_g = [zeros(1), zeros(1)]
    s_grow(outs_g, [0.6, 0.75]; idxs = [3])
    @test outs_g[1][] == only(s_grow(0.6; idxs = [3]))

    # idxs=nothing: after grow, length-3 out is accepted; short/long rejected
    out_g_ok = zeros(3)
    s_grow(out_g_ok, 0.75)
    @test out_g_ok == s_grow(0.75)
    @test_throws DimensionMismatch s_grow(zeros(2), 0.75)
    @test_throws DimensionMismatch s_grow(zeros(4), 0.75)
    outs_g_ok = [zeros(3), zeros(3)]
    s_grow(outs_g_ok, [0.6, 0.75])
    @test outs_g_ok[1] == s_grow(0.6)
    @test outs_g_ok[2] == s_grow(0.75)
    @test_throws DimensionMismatch s_grow([zeros(2), zeros(3)], [0.6, 0.75])
    @test_throws DimensionMismatch s_grow([zeros(3), zeros(4)], [0.6, 0.75])
    # before grow, length-2 out is accepted
    out_g_pre = zeros(2)
    s_grow(out_g_pre, 0.25)
    @test out_g_pre == s_grow(0.25)
    @test_throws DimensionMismatch s_grow(zeros(3), 0.25)

    # state shrunk 3 -> 2 at t = 0.5
    s_shrink = solve(
        ODEProblem(f2!, [1.0, 2.0, 3.0], (0.0, 1.0)), Tsit5();
        callback = DiscreteCallback((u, t, i) -> t == 0.5, i -> resize!(i, 2)),
        tstops = [0.5]
    )
    @test length(s_shrink(0.75)) == 2

    @test_throws BoundsError s_shrink(zeros(1), 0.75; idxs = [3])
    out_s = zeros(1)
    s_shrink(out_s, 0.25; idxs = [3])
    @test out_s[] == only(s_shrink(0.25; idxs = [3]))
    @test_throws BoundsError s_shrink(zeros(1), 0.5; idxs = [3], continuity = :right)
    out_sl = zeros(1)
    s_shrink(out_sl, 0.5; idxs = [3], continuity = :left)
    @test out_sl[] == only(s_shrink(0.5; idxs = [3], continuity = :left))
    @test_throws BoundsError s_shrink([zeros(1), zeros(1)], [0.25, 0.75]; idxs = [3])
    outs_s = [zeros(1)]
    s_shrink(outs_s, [0.25]; idxs = [3])
    @test outs_s[1][] == only(s_shrink(0.25; idxs = [3]))

    # idxs=nothing: after shrink, length-2 out is accepted; short/long rejected
    out_s_ok = zeros(2)
    s_shrink(out_s_ok, 0.75)
    @test out_s_ok == s_shrink(0.75)
    @test_throws DimensionMismatch s_shrink(zeros(1), 0.75)
    @test_throws DimensionMismatch s_shrink(zeros(3), 0.75)
    outs_s_ok = [zeros(2), zeros(2)]
    s_shrink(outs_s_ok, [0.6, 0.75])
    @test outs_s_ok[1] == s_shrink(0.6)
    @test outs_s_ok[2] == s_shrink(0.75)
    @test_throws DimensionMismatch s_shrink([zeros(1), zeros(2)], [0.6, 0.75])
    @test_throws DimensionMismatch s_shrink([zeros(2), zeros(3)], [0.6, 0.75])
    # before shrink, length-3 out is accepted
    out_s_pre = zeros(3)
    s_shrink(out_s_pre, 0.25)
    @test out_s_pre == s_shrink(0.25)
    @test_throws DimensionMismatch s_shrink(zeros(2), 0.25)
end

@testset "uninitialized scalar batch buffers" begin
    # Replace-slot path must not read vals[j] before `_set_val!` (UndefRefError)
    for initial in (1.0, big"1.0")
        sol = solve(ODEProblem((u, p, t) -> u, initial, (0.0, 1.0)), Tsit5())
        for T in (typeof(initial), Any)
            @testset "$T from $(typeof(initial))" begin
                out = Vector{T}(undef, 2)
                sol(out, [0.3, 0.7])
                @test out ≈ [sol(0.3), sol(0.7)]
            end
        end
    end
end
