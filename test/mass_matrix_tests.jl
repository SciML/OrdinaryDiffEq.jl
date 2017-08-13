using OrdinaryDiffEq, Base.Test

A = [-2.0 1 4
      4 -2 1
      2 1 3]

function f(t,u,du)
      A_mul_B!(du,A,u)
end
function f(::Type{Val{:analytic}},t,u0)
      u0.*exp(t)
end

function g(t,u,du)
      du.=u
end
function g(::Type{Val{:analytic}},t,u0)
      u0.*exp(t)
end
prob2 = ODEProblem(g,ones(3),(0.0,1.0))
prob = ODEProblem(f,ones(3),(0.0,1.0),mass_matrix=A)

######################################### Test each method for exactness

sol = solve(prob, Rosenbrock23())
sol2 = solve(prob2,Rosenbrock23())

@test norm(sol .- sol2) ≈ 0 atol=1e-11

sol = solve(prob, Rosenbrock32())
sol2 = solve(prob2,Rosenbrock32())

@test norm(sol .- sol2) ≈ 0 atol=1e-11
