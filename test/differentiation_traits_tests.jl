using OrdinaryDiffEq, Test

jac_called = Ref(false)
tgrad_called = Ref(false)

function Lotka(du,u,p,t)
  du[1] = u[1] - u[1] * u[2] # REPL[7], line 3:
  du[2] = -3 * u[2] + 1 * u[1] * u[2]
  nothing
end

function Lotka_jac(J,u,p,t)
  jac_called.x = true
  J[1,1] = 1.0 - u[2]
  J[1,2] = -u[1]
  J[2,1] = 1 * u[2]
  J[2,2] = -3 + u[1]
  nothing
end

function Lotka_tgrad(grad,u,p,t)
  tgrad_called.x = true
  grad[1] = 1 * 0
  grad[2] = 1 * 0
  nothing
end

Lotka_f = ODEFunction(Lotka,jac=Lotka_jac,tgrad=Lotka_tgrad)
prob = ODEProblem(Lotka_f,ones(2),(0.0,10.0))

good_sol = solve(prob,Rosenbrock23())

@test jac_called[]
@test tgrad_called[]

prob2 = ODEProblem(Lotka,ones(2),(0.0,10.0))

sol = solve(prob2,Rosenbrock23(autodiff=true))
@test ≈(good_sol[end],sol[end],rtol=1e-2)

sol = solve(prob2,Rosenbrock23(autodiff=true,chunk_size=1))
@test ≈(good_sol[end],sol[end],rtol=1e-2)

sol = solve(prob2,Rosenbrock23(autodiff=false))
@test ≈(good_sol[end],sol[end],rtol=1e-2)

invW_called = Ref(false)

function Lotka_invW(iW,u,p,γ,t)
  invW_called.x = true
  iW[1,1] = (1 - (1 * u[1] * u[2] * γ ^ 2) / (((1 - (-3 + u[1]) * γ) + (1 * u[1] * u[2] * γ ^ 2) / (1 - (1 - 1 * u[2]) * γ)) * (1 - (1 - 1 * u[2]) * γ))) / (1 - (1 - 1 * u[2]) * γ)
  iW[1,2] = (-(1) * u[1] * γ) / (((1 - (-3 + u[1]) * γ) + (1 * u[1] * u[2] * γ ^ 2) / (1 - (1 - 1 * u[2]) * γ)) * (1 - (1 - 1 * u[2]) * γ))
  iW[2,1] = (u[2] * γ) / (((1 - (-3 + u[1]) * γ) + (1 * u[1] * u[2] * γ ^ 2) / (1 - (1 - 1 * u[2]) * γ)) * (1 - (1 - 1 * u[2]) * γ))
  iW[2,2] = ((1 - (-3 + u[1]) * γ) + (1 * u[1] * u[2] * γ ^ 2) / (1 - (1 - 1 * u[2]) * γ)) ^ -1
  nothing
end

Lotka_f2 = ODEFunction(Lotka,jac=Lotka_jac,tgrad=Lotka_tgrad,invW = Lotka_invW)
prob3 = ODEProblem(Lotka_f2,ones(2),(0.0,10.0))
inv_sol = solve(prob3,Rosenbrock23())

@test invW_called[]
@test ≈(good_sol[end],inv_sol[end],rtol=1e-2)
