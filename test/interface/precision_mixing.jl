using OrdinaryDiffEq
function ODE(du, u, t, R, K)
    du .= u
end
params = BigFloat[1.0 0.91758707304098; 1.48439909482661 1.0]
u0 = BigFloat[0.1, 0.1]
tspan = (1.0, 31.0)
R = BigFloat[0.443280390004304303, 1.172917082211452]
K = BigFloat[13.470600276901400604, 12.52980757005]
ODE_ = (du, u, params, t) -> ODE(du, u, t, R, K)
odeProblem = ODEProblem(ODE_, u0, tspan, params)
for alg in [AutoVern8(Rodas5(), nonstifftol = 11 / 10)
            FBDF()
            QNDF()
            Tsit5()
            Rodas5P()
            TRBDF2()
            KenCarp4()
            RadauIIA5()]
    Solution = solve(odeProblem, alg, saveat = 1, abstol = 1.e-12, reltol = 1.e-6)
end
