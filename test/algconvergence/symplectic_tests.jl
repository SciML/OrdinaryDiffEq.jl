module SymplecticTests

using Test
using OrdinaryDiffEq, DiffEqBase

const ALGOS = (
    SymplecticEuler,
    VelocityVerlet,
    VerletLeapfrog,
    PseudoVerletLeapfrog,
    McAte2,
    Ruth3,
    McAte3,
    CandyRoz4,
    McAte4,
    CalvoSanz4,
    McAte42,
    McAte5,
    Yoshida6,
    KahanLi6,
    McAte8,
    KahanLi8,
    SofSpa10,
)

function dp(p, q, pa, t)
    0p .+ pa[2]
end

function dq(p, q, pa, t)
    p .* pa[1]
end

dp(res, p, q, pa, t) = (res .= dp(p, q, pa, t))
dq(res, p, q, pa, t) = (res .= dq(p, q, pa, t))

dynode(iip, dp, dq) = DynamicalODEFunction{iip}(dp, dq)

const PARAMS = ((0.1, 0.0, 1.0, 0.0), (0.1, 1.0, 1.0, -1.0))
const IIPS = (true, false)
const TSPAN = (0.0, 1.0)

solution(t, w) = (w[2] * t + w[3], (w[2] / 2 * t + w[3]) * w[1] * t + w[4])
apa(iip::Bool, x) = iip ? vcat.(x) : x

@testset "symplectic $alg-$iip-$pa" for alg in ALGOS, iip in IIPS, pa = PARAMS

    tspan = TSPAN
    t0, t1 = tspan
    dynfun = dynode(iip, dp, dq)
    p0, q0 = apa(iip, solution(t0, pa))
    prob = DynamicalODEProblem(dynfun, p0, q0, tspan, pa)
    sol = solve(prob, alg(); dt = 0.01)
    calc = sol(t1)

    @test calc[1] ≈ solution(t1, pa)[1] rtol = 1e-3
    @test calc[2] ≈ solution(t1, pa)[2] rtol = 1e-3
end
end

