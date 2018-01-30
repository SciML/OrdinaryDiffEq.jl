using OrdinaryDiffEq, Base.Test, ParameterizedFunctions

d_alembert = @ode_def DAlembert begin
    dx = a - b*x + c*t
end a b c
function (::DAlembert)(::Type{Val{:analytic}},u0,p,t::Number)
    a,b,c = p
    ebt = exp(b*t)
    exp(-b*t)*(-a*b + c + ebt*(a*b + c*(b*t - 1)) + b^2 * u0)/(b^2)
end

p = (1., 2., 3.)
u0 = [1.0]
tspan = (0.0,10.0)
prob = ODEProblem(d_alembert,u0,tspan,p)

sol = solve(prob,Tsit5(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,Rosenbrock23(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,Rodas42(),abstol=1e-10,reltol=1e-10)
@test_broken sol.errors[:l2] < 1e-7
sol = solve(prob,Rodas5(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,TRBDF2(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,Trapezoid(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,KenCarp3(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,KenCarp4(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7

prob2 = ODEProblem((du,u,p,t)->d_alembert(du,u,p,t),u0,tspan,p)

sol = solve(prob2,Rosenbrock23(),abstol=1e-10,reltol=1e-10)
sol = solve(prob2,Rodas4(),abstol=1e-10,reltol=1e-10)
