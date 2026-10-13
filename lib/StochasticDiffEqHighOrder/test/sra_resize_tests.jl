using StochasticDiffEqHighOrder
using SciMLBase
using Test

function _step_after_resize!(integ, n = 10)
    for _ in 1:n
        step!(integ)
    end
    return integ
end

function _resize_and_step!(prob, alg; kwargs...)
    integ = init(prob, alg; dt = 0.01, kwargs...)
    step!(integ)
    resize!(integ, 3)
    integ.u[3] = 1.0
    c = integ.cache
    @test all(x -> length(x) == 3, c.H0)
    if hasfield(typeof(c), :H1)
        @test all(x -> length(x) == 3, c.H1)
    end
    _step_after_resize!(integ)
    @test length(integ.u) == 3
    @test all(isfinite, integ.u)

    deleteat!(integ, 3)
    @test all(x -> length(x) == 2, c.H0)
    step!(integ)
    @test length(integ.u) == 2

    addat!(integ, 3:3)
    integ.u[3] = 1.0
    @test all(x -> length(x) == 3, c.H0)
    step!(integ)
    @test length(integ.u) == 3
    @test all(isfinite, integ.u)
    return nothing
end

@testset "SRA/SRI resize! stage buffers (issue 4720)" begin
    f!(du, u, p, t) = (du .= -u)
    g!(du, u, p, t) = (du .= 0.1)
    prob = SDEProblem(f!, g!, ones(2), (0.0, 1.0))
    for alg in (SRA(), SRI())
        _resize_and_step!(prob, alg; adaptive = false)
    end
    _resize_and_step!(prob, SRA())
end
