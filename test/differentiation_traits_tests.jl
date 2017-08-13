using OrdinaryDiffEq, Base.Test

jac_called = false
tgrad_called = false

function Lotka(t,u,du)
  du[1] = u[1] - u[1] * u[2] # REPL[7], line 3:
  du[2] = -3 * u[2] + 1 * u[1] * u[2]
  nothing
end

function Lotka(::Type{Val{:jac}},t,u,J)
  global jac_called
  jac_called = true
  J[1,1] = 1.0 - u[2]
  J[1,2] = -u[1]
  J[2,1] = 1 * u[2]
  J[2,2] = -3 + u[1]
  nothing
end

function Lotka(::Type{Val{:tgrad}},t,u,grad)
  global tgrad_called
  tgrad_called = true
  grad[1] = 1 * 0
  grad[2] = 1 * 0
  nothing
end

prob = ODEProblem(Lotka,ones(2),(0.0,10.0))

good_sol = solve(prob,Rosenbrock23())

@test jac_called == true
@test tgrad_called == true

function Lotka2(t,u,du)
  du[1] = u[1] - u[1] * u[2] # REPL[7], line 3:
  du[2] = -3 * u[2] + 1 * u[1] * u[2]
  nothing
end

prob2 = ODEProblem(Lotka2,ones(2),(0.0,10.0))

sol = solve(prob2,Rosenbrock23(autodiff=true))
@test ≈(good_sol[end],sol[end],rtol=1e-2)

sol = solve(prob2,Rosenbrock23(autodiff=true,chunk_size=1))
@test ≈(good_sol[end],sol[end],rtol=1e-2)

sol = solve(prob2,Rosenbrock23(autodiff=false))
@test ≈(good_sol[end],sol[end],rtol=1e-2)

invW_called = false

function Lotka3(t,u,du)
  du[1] = u[1] - u[1] * u[2] # REPL[7], line 3:
  du[2] = -3 * u[2] + 1 * u[1] * u[2]
  nothing
end

function Lotka3(::Type{Val{:invW}},t,u,γ,iW)
  global invW_called
  invW_called = true
  iW[1,1] = (1 - (1 * u[1] * u[2] * γ ^ 2) / (((1 - (-3 + u[1]) * γ) + (1 * u[1] * u[2] * γ ^ 2) / (1 - (1 - 1 * u[2]) * γ)) * (1 - (1 - 1 * u[2]) * γ))) / (1 - (1 - 1 * u[2]) * γ)
  iW[1,2] = (-(1) * u[1] * γ) / (((1 - (-3 + u[1]) * γ) + (1 * u[1] * u[2] * γ ^ 2) / (1 - (1 - 1 * u[2]) * γ)) * (1 - (1 - 1 * u[2]) * γ))
  iW[2,1] = (u[2] * γ) / (((1 - (-3 + u[1]) * γ) + (1 * u[1] * u[2] * γ ^ 2) / (1 - (1 - 1 * u[2]) * γ)) * (1 - (1 - 1 * u[2]) * γ))
  iW[2,2] = ((1 - (-3 + u[1]) * γ) + (1 * u[1] * u[2] * γ ^ 2) / (1 - (1 - 1 * u[2]) * γ)) ^ -1
  nothing
end

function Lotka3(::Type{Val{:tgrad}},t,u,grad)
  global tgrad_called
  tgrad_called = true
  grad[1] = 1 * 0
  grad[2] = 1 * 0
  nothing
end

prob3 = ODEProblem(Lotka3,ones(2),(0.0,10.0))
inv_sol = solve(prob3,Rosenbrock23())

@test invW_called
@test ≈(good_sol[end],inv_sol[end],rtol=1e-2)
