using OrdinaryDiffEq, Plots
srand(100)
setprecision(400)

f = (t,u) -> (2u)
analytic = (t,u0) -> u0*exp(t)
prob = ODETestProblem(f,1/2,analytic)


sol3 =solve(prob,RK4(),dt=1/2^(6),abstol=1,reltol=0)

prob = ODETestProblem(f,BigInt(1)//BigInt(2),analytic)

sol =solve(prob,RK4(),dt=BigInt(1)//BigInt(2)^(6),abstol=1,reltol=0)
sol2 =solve(prob,RK4(),dt=BigInt(1)/BigInt(2)^(6),abstol=1,reltol=0)

sol.u
sol2.u
sol3.u
@test 6.37e-16 < abs(sol.u[end] - sol3.u[end]) < 6.38e-16
@test 6.37e-16  < abs(sol2.u[end] - sol3.u[end]) <6.38e-16
@test sol2.u[end] - sol.u[end] < big(1.73e-77)
@test eltype(sol.u) == Rational{BigInt}
@test eltype(sol2.u) == Rational{BigInt}
@test eltype(sol3.u) == Float64

sol4 =solve(prob,DP5(),dt=BigInt(1)//BigInt(2)^(3),adaptive=false)

@test eltype(sol4.u) == Rational{BigInt}

tab = constructDormandPrince8_64bit(Rational{BigInt})
sol5 =solve(prob,ExplicitRK(),dt=BigInt(1)//BigInt(2)^(3),abstol=1,reltol=0,tableau=tab,adaptive=false)

@test 5.72e-8 < abs(float(sol5.u[end]) - sol3.u[end]) < 5.73e-8
@test typeof(sol5.u[end]) == Rational{BigInt}

prob = ODETestProblem(f,1/2+1/2im,analytic)
