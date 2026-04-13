using Test, LinearAlgebra
using OrdinaryDiffEqSymplecticRK, DiffEqBase
using OrdinaryDiffEqRKN

# algorithm, dq(p) != p, convergence order
const ALGOS = (
    (SymplecticEuler, true, 1),
    (VelocityVerlet, false, 2),
    (VerletLeapfrog, true, 2),
    (LeapfrogDriftKickDrift, true, 2),
    (PseudoVerletLeapfrog, true, 2),
    (McAte2, true, 2),
    (Ruth3, true, 3),
    (McAte3, true, 3),
    (CandyRoz4, true, 4),
    (McAte4, true, 4),
    (CalvoSanz4, true, 4),
    (McAte42, true, 1), # known to be broken
    (McAte5, true, 5),
    (Yoshida6, true, 6),
    (KahanLi6, true, 6),
    (McAte8, true, 8),
    (KahanLi8, true, 8),
    (SofSpa10, true, 10),
)

function dp(p, q, pa, t)
    return 0q .+ pa[2]
end

function dq(p, q, pa, t)
    return p .* pa[1]
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
    return println()
end

@testset "symplectic $alg-$iip-$pa" for (alg, x, d) in ALGOS, iip in IIPS, pa in PARAMS
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

function motionfuncDirect1(dv, v, u, p, t)
    # 1:Electron, 2: Be
    ω_1, ω_2, γ, m_1, m_2, η, ω_d = p
    dv[1] = -ω_1^2 * u[1] * (1 + η * cos(ω_d * t)) - γ * u[2] / m_1
    return dv[2] = -ω_2^2 * u[2] - γ * u[1] / m_2
end

function motionfuncDirect1(v, u, p, t)
    # 1:Electron, 2: Be
    ω_1, ω_2, γ, m_1, m_2, η, ω_d = p
    return [
        -ω_1^2 * u[1] * (1 + η * cos(ω_d * t)) - γ * u[2] / m_1,
        -ω_2^2 * u[2] - γ * u[1] / m_2,
    ]
end

param = [
    90386.15717208837, 3938.9288690708827, 8560.718748264337, 0.000544617021484666,
    8.947079933513658, 0.7596480420227258, 78778.57738141765,
]
u0_direct = zeros(2) # mm, mm
v0_direct = [0.0, 135.83668926684385]
tspan = (0.0, 1.321179076090661)
prob_direct = SecondOrderODEProblem(motionfuncDirect1, v0_direct, u0_direct, tspan, param)
dt = 2.0e-8
ref = solve(
    prob_direct, DPRKN12(), abstol = 1.0e-12, reltol = 1.0e-12, maxiters = 1.0e7, saveat = 0.01
)

@testset "symplectic time-dependent $alg" for (alg, x, d) in ALGOS
    sol = solve(prob_direct, alg(), dt = dt, saveat = 0.01)
    if alg <: Yoshida6
        @test maximum(ref[4, :] - sol[4, :]) < 9.0e-3
    else
        @test maximum(ref[4, :] - sol[4, :]) < 3.0e-3
    end
end
