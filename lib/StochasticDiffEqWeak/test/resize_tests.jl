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
