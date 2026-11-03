using StochasticDiffEqHighOrder, Test
using SciMLBase: full_cache, addat!

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

function seeded_resize_run(alg)
    f(du, u, p, t) = (du .= -0.5 .* u)
    g(du, u, p, t) = (du .= 0.1 .* u)
    prob = SDEProblem(f, g, ones(2), (0.0, 1.0))
    integ = init(prob, alg; dt = 0.01, adaptive = false, seed = 99)
    foreach(_ -> step!(integ), 1:3)
    resize!(integ, 5)
    integ.u[3:5] .= 1.0
    foreach(_ -> step!(integ), 1:3)
    resize!(integ, 3)
    foreach(_ -> step!(integ), 1:2)
    deleteat!(integ, 1)
    step!(integ)
    addat!(integ, 2:2)
    integ.u[2] = 1.0
    foreach(_ -> step!(integ), 1:2)
    return integ.u
end

# Reference values pin the order of the noise draws made when the state grows.
@testset "$(nameof(typeof(alg))) seeded trajectory through resize!/deleteat!/addat!" for (alg, ref) in (
        (SRIW1(), [0.9365594868584296, 0.9951144286306745, 0.9993877804414055]),
        (SRA1(), [0.9365645455676647, 0.9950617203015882, 0.9993902725961651]),
    )
    @test seeded_resize_run(alg) ≈ ref rtol = 1.0e-12
end
