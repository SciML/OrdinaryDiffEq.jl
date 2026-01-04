using OrdinaryDiffEq, Test

function f(du, u, p, t)
    Gut, Cent, Periph, Resp = u
    Ka, CL, Vc, Q, Vp, Kin, Kout, IC50, IMAX, γ = p
    dGut = -Ka * Gut
    dCent = Ka * Gut - CL * (Cent / Vc) - Q * (Cent / Vc) + Q * (Periph / Vp)
    dPeriph = Q * (Cent / Vc) - Q * (Periph / Vp)
    dResp = Kin * (1 - (IMAX * (Cent / Vc)^γ / (IC50^γ + (Cent / Vc)^γ))) - Kout * Resp
    du[1] = dGut
    du[2] = dCent
    du[3] = dPeriph
    return du[4] = dResp
end

p = (
    Ka = 2.769674683505292, CL = 0.8315548845614537, Vc = 45.87339297906569,
    Q = 1.4544520084387496, Vp = 9.460961128900388, Kin = 6.360418843010778,
    Kout = 3.524234094210465, IC50 = 0.4483840673123028, IMAX = 1.0, γ = 2.0,
)

u0 = [
    0.35793622925128044,
    0.7123008259214687,
    0.10685666303615937,
    0.0,
]
prob = ODEProblem(f, u0, (0, 26.0), p)
sol = solve(prob, AutoTsit5(Rodas5(), dtfac = 2))
@test sol([21.2])[4, 1] ≈ 1.803 atol = 1.0e-2

#=
using Plots
plot(sol, vars=4, lab="dtfac=1")
sol = solve(prob, AutoTsit5(Rodas5(), dtfac=2))
plt1 = plot!(sol, vars=4, lab="dtfac=2")
plt2 = plot(sol.t, sol.alg_choice, lab="alg_choice")
plot(plt1, plt2, layout=(2,1))
plot!(;title="AutoTsit5(Rodas5())",dpi=200)
=#
