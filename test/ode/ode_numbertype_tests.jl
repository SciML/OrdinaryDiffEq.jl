using OrdinaryDiffEq, Plots
srand(100)
setprecision(400)

f = (t,u) -> (2u)
analytic = (t,u0) -> u0*exp(t)
prob = ODETestProblem(f,1/2,analytic)


sol3 =solve(prob,RK4,dt=1/2^(6),abstol=1,reltol=0)

prob = ODETestProblem(f,BigInt(1)//BigInt(2),analytic)

sol =solve(prob,RK4,dt=BigInt(1)//BigInt(2)^(6),abstol=1,reltol=0)
sol2 =solve(prob,RK4,dt=BigInt(1)/BigInt(2)^(6),abstol=1,reltol=0)

sol.u
sol2.u
sol3.u
bool1 = 6.37e-16 < abs(sol.u[end] - sol3.u[end]) < 6.38e-16
bool2 = 6.37e-16  < abs(sol2.u[end] - sol3.u[end]) <6.38e-16
bool3 = sol2.u[end] - sol.u[end] < big(1.73e-77)
bool4 = eltype(sol.u) == Rational{BigInt}
bool5 = eltype(sol2.u) == Rational{BigInt}
bool6 = eltype(sol3.u) == Float64

sol4 =solve(prob,DP5,dt=BigInt(1)//BigInt(2)^(3),adaptive=false)

bool9 = eltype(sol4.u) == Rational{BigInt}

tab = constructDormandPrince8_64bit(Rational{BigInt})
sol5 =solve(prob,ExplicitRK,dt=BigInt(1)//BigInt(2)^(3),abstol=1,reltol=0,tableau=tab,adaptive=false)

bool7 = 5.72e-8 < abs(float(sol5.u[end]) - sol3.u[end]) < 5.73e-8
bool8 = typeof(sol5.u[end]) == Rational{BigInt}

prob = ODETestProblem(f,1/2+1/2im,analytic)

bool1 && bool2 && bool3 && bool4 && bool5 && bool6 && bool7 && bool8 && bool9
