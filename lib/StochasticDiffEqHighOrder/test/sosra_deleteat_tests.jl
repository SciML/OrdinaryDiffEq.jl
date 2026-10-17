using StochasticDiffEqHighOrder, Test
using SciMLBase: full_cache

f(du, u, p, t) = (du .= -u)
g(du, u, p, t) = (du .= 0.1)
prob = SDEProblem(f, g, ones(2), (0.0, 1.0))

@testset "$(nameof(typeof(alg))) resize! then deleteat!" for alg in (SRA2(), SRA3(), SOSRA(), SOSRA2())
    integ = init(prob, alg)
    @test integ.cache.tmp !== integ.cache.k1
    resize!(integ, 3)
    @test all(length(c) == 3 for c in full_cache(integ.cache))
    step!(integ)
    @test length(integ.u) == 3
    deleteat!(integ, 3)
    @test all(length(c) == 2 for c in full_cache(integ.cache))
    step!(integ)
    @test length(integ.u) == 2
    @test all(isfinite, integ.u)
end
