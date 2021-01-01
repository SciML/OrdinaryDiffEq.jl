
using Test, LinearAlgebra
using OrdinaryDiffEq, DiffEqBase

const ALGOS = (
    (SymplecticEuler, true, 1),
    (VelocityVerlet, false, 2),
    (VerletLeapfrog, true, 2),
    (PseudoVerletLeapfrog, true, 2),
    (McAte2, true, 2),
    (Ruth3, true, 3),
    (McAte3, true, 3),
    (CandyRoz4, true, 4),
    (McAte4, true, 4),
    (CalvoSanz4, true, 4),
    # (McAte42, true, 4), known to be broken
    (McAte5, true, 5),
    (Yoshida6, true, 6),
    (KahanLi6, true, 6),
    (McAte8, true, 8),
    (KahanLi8, true, 8),
    (SofSpa10, true, 10),
)

function dp(p, q, pa, t)
    0q .+ pa[2]
end

function dq(p, q, pa, t)
    p .* pa[1]
end

dp(res, p, q, pa, t) = (res .= dp(p, q, pa, t))
dq(res, p, q, pa, t) = (res .= dq(p, q, pa, t))

dynode(iip, dp, dq) = DynamicalODEFunction{iip}(dp, dq)

# [0:1] used in dp, dq; [3:4] start values for p0, q0
const PARAMS = ((1.0, 0.1, 1.0, 0.0), (0.1, 1.0, 1.0, -1.0))
const IIPS = (true, false)
const TSPAN = (0.0, 1.0)

solution(t, w) = (w[2] * t + w[3], (w[2] / 2 * t + w[3]) * w[1] * t + w[4])
apa(iip::Bool, x) = iip ? vcat.(x) : x
errorbound(dt, d, x) = 100 * abs(dt)^d + 1000 * eps(norm(x))
function printerrors(text, calc, solution, pa, t1)
    print(text, ": ")
    print(norm(calc[1] - solution(t1, pa)[1]), " ")
    print(norm(calc[2] - solution(t1, pa)[2]))
    println()
end

@testset "symplectic $alg-$iip-$pa" for (alg, x, d) in ALGOS, iip in IIPS, pa = PARAMS

    dt = 0.01
    tspan = TSPAN
    t0, t1 = tspan
    dynfun = dynode(iip, dp, dq)
    p0, q0 = apa(iip, solution(t0, pa))
    prob = DynamicalODEProblem(dynfun, p0, q0, tspan, pa)

    if x || pa[1] == 1
        sol = solve(prob, alg(); dt = dt)    
        calc = sol(t1)
        # printerrors("$alg-$iip-$pa", calc, solution, pa, t1)
        @test calc[1] ≈ solution(t1, pa)[1] rtol = errorbound(dt, d, calc[1])
        @test calc[2] ≈ solution(t1, pa)[2] rtol = errorbound(dt, d, calc[2])
    else
        @test_throws ArgumentError solve(prob, alg(); dt = dt)
    end
end

