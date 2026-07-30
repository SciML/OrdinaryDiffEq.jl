using StochasticDiffEqROCK, DiffEqDevTools, Test
using LinearAlgebra, Random

# SROCK2 is weak order 2 for general (non-diagonal) Ito SDEs, see Abdulle,
# Vilmart & Zygalakis, SIAM J. Sci. Comput. 35(4), 2013, Theorem 3.4. The m = 1
# problems guard the #3170 routing fix and the m = 2 problem covers the
# general-noise branch of perform_step!. Every B_k of both problems is diagonal,
# so their diffusions commute and the J_{q,r} Rademacher cross terms of eq. (17)
# cancel: reverting the eq. (17) hunks leaves all four order estimates below
# unchanged. The non-commutative second moment check at the bottom is what
# covers those terms.

const drift_c = -0.5

# n = 2 states driven by the same single Brownian motion (m = 1).
const diff_c = 0.1
f1_oop(u, p, t) = drift_c .* u
g1_oop(u, p, t) = diff_c .* reshape(u, 2, 1)
f1_iip(du, u, p, t) = (du .= drift_c .* u)
g1_iip(du, u, p, t) = (du .= diff_c .* u)
function analytic1(u0, p, t, W)
    return u0 .* exp.((drift_c - diff_c^2 / 2) * t .+ diff_c .* W[1])
end

prob1_oop = SDEProblem(
    SDEFunction(f1_oop, g1_oop; analytic = analytic1),
    [1.0, 1.0], (0.0, 1.0); noise_rate_prototype = zeros(2, 1)
)
prob1_iip = SDEProblem(
    SDEFunction(f1_iip, g1_iip; analytic = analytic1),
    [1.0, 1.0], (0.0, 1.0); noise_rate_prototype = zeros(2, 1)
)

# n = 2 states, each driven by both of m = 2 Brownian motions.
const diff_B = [0.2 0.1; 0.1 0.2]
f2_oop(u, p, t) = drift_c .* u
g2_oop(u, p, t) = diff_B .* u
f2_iip(du, u, p, t) = (du .= drift_c .* u)
g2_iip(du, u, p, t) = (du .= diff_B .* u)
function analytic2(u0, p, t, W)
    drift = drift_c .- vec(sum(abs2, diff_B, dims = 2)) ./ 2
    return u0 .* exp.(drift .* t .+ diff_B * W[1:2])
end

prob2_oop = SDEProblem(
    SDEFunction(f2_oop, g2_oop; analytic = analytic2),
    [1.0, 1.0], (0.0, 1.0); noise_rate_prototype = zeros(2, 2)
)
prob2_iip = SDEProblem(
    SDEFunction(f2_iip, g2_iip; analytic = analytic2),
    [1.0, 1.0], (0.0, 1.0); noise_rate_prototype = zeros(2, 2)
)

dts = 1 .// 2 .^ (6:-1:2)

Random.seed!(100)
println("SROCK2 non-diagonal m=1")
for prob in (prob1_oop, prob1_iip)
    @time sim = test_convergence(
        dts, prob, SROCK2(), trajectories = Int(5.0e4),
        save_everystep = false, weak_timeseries_errors = false
    )
    @show sim.𝒪est[:weak_final]
    @test abs(sim.𝒪est[:weak_final] - 2) < 0.5
end

println("SROCK2 non-diagonal m=2")
for prob in (prob2_oop, prob2_iip)
    @time sim = test_convergence(
        dts, prob, SROCK2(), trajectories = Int(5.0e4),
        save_everystep = false, weak_timeseries_errors = false
    )
    @show sim.𝒪est[:weak_final]
    @test abs(sim.𝒪est[:weak_final] - 2) < 0.5
end

# Non-commutative diffusion g(u) = [0 a*u2; b*u1 0], i.e. the two noise channels
# are B1 = [0 0; b 0] and B2 = [0 a; 0 0] with B1*B2 != B2*B1, so the eq. (17)
# cross terms no longer cancel. For a linear Ito SDE the second moment
# P(t) = E[u u'] solves P' = A*P + P*A' + sum_k B_k*P*B_k', so that
# vec(P)(T) = exp(M*T)*vec(u0*u0') with M = kron(I, A) + kron(A, I) +
# sum_k kron(B_k, B_k). The cross entry E[u1*u2](T) is the one the cross terms
# move: at dt = 1/2 with 2e5 trajectories its error stays under 9e-3 with
# eq. (17) active and sits near 6e-2 with those terms commented out, so the
# tolerance below is a regression guard for them rather than an accuracy claim.
const a_nc = 1.2
const b_nc = 1.0
const A_nc = [-0.5 0.0; 0.0 -0.5]
const B1_nc = [0.0 0.0; b_nc 0.0]
const B2_nc = [0.0 a_nc; 0.0 0.0]
@test B1_nc * B2_nc != B2_nc * B1_nc

fnc_oop(u, p, t) = A_nc * u
gnc_oop(u, p, t) = [0.0 a_nc * u[2]; b_nc * u[1] 0.0]
fnc_iip(du, u, p, t) = mul!(du, A_nc, u)
function gnc_iip(du, u, p, t)
    du[1, 1] = 0.0
    du[1, 2] = a_nc * u[2]
    du[2, 1] = b_nc * u[1]
    du[2, 2] = 0.0
    return nothing
end

u0_nc = [1.0, 1.0]
T_nc = 1.0
M_nc = let Id = Matrix(1.0I, 2, 2)
    kron(Id, A_nc) + kron(A_nc, Id) + kron(B1_nc, B1_nc) + kron(B2_nc, B2_nc)
end
P_nc = reshape(exp(M_nc * T_nc) * vec(u0_nc * u0_nc'), 2, 2)

prob_nc_oop = SDEProblem(
    fnc_oop, gnc_oop, u0_nc, (0.0, T_nc); noise_rate_prototype = zeros(2, 2)
)
prob_nc_iip = SDEProblem(
    fnc_iip, gnc_iip, u0_nc, (0.0, T_nc); noise_rate_prototype = zeros(2, 2)
)

println("SROCK2 non-commutative m=2 second moment")
for prob in (prob_nc_oop, prob_nc_iip)
    Random.seed!(100)
    N_nc = 200000
    P = zeros(2, 2)
    @time for _ in 1:N_nc
        u = solve(
            prob, SROCK2(); dt = 1 / 2, adaptive = false,
            save_everystep = false, save_start = false
        ).u[end]
        P .+= u * u'
    end
    P ./= N_nc
    @show abs(P[1, 2] - P_nc[1, 2])
    @test abs(P[1, 2] - P_nc[1, 2]) < 2.0e-2
end
