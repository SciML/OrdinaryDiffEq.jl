using DelayDiffEq, Test
using OrdinaryDiffEqVerner
using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEqHighOrderRK
using OrdinaryDiffEqFeagin
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using Random
Random.seed!(213)

const CACHE_TEST_ALGS = [
    Vern7(), Euler(), Midpoint(), RK4(), SSPRK22(), SSPRK33(),
    ORK256(), CarpenterKennedy2N54(), SHLDDRK64(), DGLDDRK73_C(),
    CFRLDDRK64(), TSLDDRK74(),
    CKLLSRK43_2(),
    ParsaniKetchesonDeconinck3S32(),
    BS3(), BS5(), DP5(), DP8(), Feagin10(), Feagin12(), Feagin14(), TanYam7(),
    Tsit5(), TsitPap8(), Vern6(), Vern7(), Vern8(), Vern9(), OwrenZen3(), OwrenZen4(),
    OwrenZen5(),
]

function f(du, u, h, p, t)
    for i in 1:length(u)
        du[i] = (0.3 / length(u)) * u[i]
    end
    return nothing
end

condition(u, t, integrator) = 1 - maximum(u)

function affect!(integrator)
    u = integrator.u
    resize!(integrator, length(u) + 1)
    maxidx = findmax(u)[2]
    Θ = rand() / 5 + 0.25
    u[maxidx] = Θ
    u[end] = 1 - Θ
    return nothing
end

const prob = DDEProblem(
    f, [0.2], nothing, (0.0, 10.0);
    callback = ContinuousCallback(condition, affect!)
)

println("Check for stochastic errors")
for i in 1:10
    @test_nowarn solve(prob, MethodOfSteps(Tsit5()); dt = 0.5)
end

println("Check some other integrators")

println("Rosenbrock23")
@test_nowarn solve(prob, MethodOfSteps(Rosenbrock23()); dt = 0.5)

@testset "$(nameof(typeof(alg)))" for alg in CACHE_TEST_ALGS
    @test_nowarn solve(prob, MethodOfSteps(alg); dt = 0.5)
end

println("Rodas4")
@test_nowarn solve(prob, MethodOfSteps(Rodas4()); dt = 0.5)

println("Rodas5")
@test_nowarn solve(prob, MethodOfSteps(Rodas5()); dt = 0.5)
