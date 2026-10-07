using StochasticDiffEqWeak, Test
using SciMLBase: addat!

f(du, u, p, t) = (du .= -0.5 .* u)
g(du, u, p, t) = (du .= 0.1 .* u)
prob = SDEProblem(f, g, ones(2), (0.0, 1.0))

const STAGE_FIELDS = (:g2, :g3, :g4, :H12, :H13, :H14, :H22, :H23, :Yp, :Ym)

function check_stage_lengths(cache, n)
    for name in STAGE_FIELDS
        hasproperty(cache, name) || continue
        v = getproperty(cache, name)
        v isa AbstractVector{<:AbstractArray} || continue
        @test length(v) == n
        @test all(length.(v) .== n)
    end
    return
end

function step_and_check!(integ, n)
    check_stage_lengths(integ.cache, n)
    step!(integ)
    @test length(integ.u) == n
    return @test all(isfinite, integ.u)
end

@testset "$(nameof(typeof(alg))) resize!/deleteat!/addat! with diagonal noise" for alg in (
        DRI1(), RI1(), RDI2WM(), RS1(), RS2(), W2Ito1(), IRI1(),
    )
    integ = init(prob, alg; dt = 0.01, adaptive = false)
    step!(integ)

    resize!(integ, 4)
    integ.u[3:4] .= 1.0
    step_and_check!(integ, 4)
    @test integ.u[3] != 1.0

    resize!(integ, 3)
    step_and_check!(integ, 3)

    # IRI1's nonlinear solver buffers are only resized by `resize!`
    alg isa IRI1 && continue

    deleteat!(integ, 1)
    step_and_check!(integ, 2)

    addat!(integ, 3:3)
    integ.u[3] = 1.0
    step_and_check!(integ, 3)
end

# PL1WM's Z process holds one entry per pair of noise dimensions, W2Ito1's always two.
@testset "$(nameof(typeof(alg))) (adaptive = $adaptive) Z process follows the state" for (alg, zlen, adaptive) in (
        (PL1WM(), m -> m * (m - 1) ÷ 2, false), (W2Ito1(), m -> 2, false), (W2Ito1(), m -> 2, true),
    )
    integ = init(prob, alg; dt = 0.01, adaptive)
    step!(integ)
    for (n, change!) in (
            (4, i -> resize!(i, 4)), (3, i -> resize!(i, 3)), (5, i -> resize!(i, 5)),
            (4, i -> deleteat!(i, 2)), (5, i -> addat!(i, 1:1)), (2, i -> resize!(i, 2)),
        )
        change!(integ)
        integ.u .= 1.0
        @test length(integ.W.dZ) == zlen(n)
        @test length(integ.cache._dZ) == zlen(n)
        check_stage_lengths(integ.cache, n)
        step!(integ)
        step!(integ)
        @test length(integ.u) == n
        @test all(isfinite, integ.u)
    end
    adaptive && @test all(c -> length(c[3]) == zlen(2), integ.W.S₂.data)
end
