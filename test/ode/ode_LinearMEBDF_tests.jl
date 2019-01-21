using OrdinaryDiffEq, Test, DiffEqDevTools,  DiffEqOperators, SparseArrays, LinearAlgebra
alg = LinearMEBDF(linsolve=LinSolveFactorize(lu))
dts = 0.5 .^(12:-1:7)
testTol = 0.2
u0 = rand(100)
tspan = (0.0, 1.0)

###constant linear Opeartor:
A = spdiagm(-1 => ones(99) , 0 => fill(-2, 100) , 1 => ones(99));
f = DiffEqArrayOperator(A);

sol_analytic = (u0,p,t) -> exp(convert(Array,(t * A)))*u0
prob = ODEProblem(ODEFunction(f; analytic = sol_analytic),u0,tspan)

sim = test_convergence(dts, prob, alg)
println("LinearMEBDF, order = ", sim.𝒪est[:l2])
@test sim.𝒪est[:final] ≈ 2 atol=0.2


### M≠I:
M = rand(100,100)
sol_analytic = (u0,p,t) -> exp(convert(Array,(t*inv(M)*A))) * u0
prob = ODEProblem(ODEFunction(f;
   analytic = sol_analytic, mass_matrix = M),u0,tspan)

sim = test_convergence(dts, prob, alg)
println("LinearMEBDF, order = ", sim.𝒪est[:l2])
@test sim.𝒪est[:final] ≈ 2 atol=0.2

####Non constant linear Operator:

B = spdiagm(0 => ones(100));
update_func = (_A,u,p,t) -> _A.nzval  .= t
L = DiffEqArrayOperator(B; update_func = update_func);
sol_analytic = (u0,p,t) -> exp(t^2/2) .* u0
prob = ODEProblem(ODEFunction(L; analytic = sol_analytic),u0,tspan)

sim = test_convergence(dts, prob, alg)
println("LinearMEBDF, order = ", sim.𝒪est[:l2])
@test sim.𝒪est[:final] ≈ 2 atol=0.2

### M≠I:
M= 0.5*Diagonal(ones(100))
C = spdiagm(0 => ones(100));
update_func = (_A,u,p,t) -> _A.nzval  .= t
  L = DiffEqArrayOperator(C; update_func = update_func);

  sol_analytic = (u0,p,t) -> exp(t^2/(2*0.5)) .* u0
  prob = ODEProblem(ODEFunction(L;
     analytic = sol_analytic, mass_matrix = M),u0,tspan)

sim = test_convergence(dts, prob, alg)
println("LinearMEBDF, order = ", sim.𝒪est[:l2])
@test sim.𝒪est[:final] ≈ 2 atol=0.2
