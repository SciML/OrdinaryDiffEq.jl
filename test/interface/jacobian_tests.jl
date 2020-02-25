using OrdinaryDiffEq, Test

function d_alembert(du,u,p,t)
  du[1] = p[1] - p[2]*u[1] + p[3]*t
end

function d_alembert_jac(J,u,p,t)
  J[1] = -p[2]
end

function d_alembert_analytic(u0,p,t::Number)
    a,b,c = p
    ebt = exp(b*t)
    @. exp(-b*t)*(-a*b + c + ebt*(a*b + c*(b*t - 1)) + b^2 * u0)/(b^2)
end

p = (1., 2., 3.)
u0 = [1.0]
tspan = (0.0,10.0)
prob = ODEProblem(ODEFunction(d_alembert,
                              jac = d_alembert_jac,
                              analytic=d_alembert_analytic),
                              u0,tspan,p)

sol = solve(prob,Tsit5(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,Rosenbrock23(),abstol=1e-8,reltol=1e-8)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,Rodas4(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,Veldd4(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,Rodas5(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,TRBDF2(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 2e-6
sol = solve(prob,Trapezoid(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 2e-6
sol = solve(prob,KenCarp3(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 8e-4
sol = solve(prob,KenCarp4(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7

using ModelingToolkit
function lotka(du,u,p,t)
  x = u[1]
  y = u[2]
  du[1] = p[1]*x - p[2]*x*y
  du[2] = -p[3]*y + p[4]*x*y
end

prob = ODEProblem(lotka,[1.0,1.0],(0.0,1.0),[1.5,1.0,3.0,1.0])
de = ModelingToolkit.modelingtoolkitize(prob)
prob2 = remake(prob, f=ODEFunction(de...; jac=true, Wfact=true))

sol = solve(prob, TRBDF2())

for Alg in [Rodas5, Rosenbrock23, TRBDF2, KenCarp4]
  @test Array( solve(prob2, Alg(), tstops=sol.t, adaptive=false) ) ≈ Array( solve(prob, Alg(), tstops=sol.t, adaptive=false) ) atol=1e-4
end

_Wfact = eval(ModelingToolkit.generate_factorized_W(de[1])[1][2])
_Wfact_t = eval(ModelingToolkit.generate_factorized_W(de[1])[2][2])

lotka_only_Wfact   = remake(prob,f=ODEFunction(lotka,Wfact=_Wfact))
lotka_only_Wfact_t = remake(prob,f=ODEFunction(lotka,Wfact_t=_Wfact_t))

sol = solve(lotka_only_Wfact, TRBDF2())
sol = solve(lotka_only_Wfact_t, TRBDF2())
sol = solve(lotka_only_Wfact, Rosenbrock23())
sol = solve(lotka_only_Wfact_t, Rosenbrock23())
