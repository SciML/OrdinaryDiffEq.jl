using OrdinaryDiffEq, Base.Test

const A = [-2.0 1 4
            4 -2 1
            2 1 3]
const b = A*ones(3)
function f(t,u,du)
      A_mul_B!(du,A,u)
      tmp = t*b
      du .+= tmp
end
function f(::Type{Val{:analytic}},t,u0)
      @. 2ones(3)*exp(t) - t - 1
end
function g(t,u,du)
      du .= u + t
end
function g(::Type{Val{:analytic}},t,u0)
      @. 2ones(3)*exp(t) - t - 1
end
prob2 = ODEProblem(g,ones(3),(0.0,1.0))
prob = ODEProblem(f,ones(3),(0.0,1.0),mass_matrix=A)

######################################### Test each method for exactness

sol = solve(prob  ,Rosenbrock23())
sol2 = solve(prob2,Rosenbrock23())

@test norm(sol .- sol2) ≈ 0 atol=1e-11

sol = solve(prob, Rosenbrock32())
sol2 = solve(prob2,Rosenbrock32())

@test norm(sol .- sol2) ≈ 0 atol=1e-11

sol = solve(prob,  ROS3P())
sol2 = solve(prob2,ROS3P())

@test norm(sol .- sol2) ≈ 0 atol=1e-11

sol = solve(prob,  Rodas3())
sol2 = solve(prob2,Rodas3())

@test norm(sol .- sol2) ≈ 0 atol=1e-11

sol = solve(prob,  RosShamp4())
sol2 = solve(prob2,RosShamp4())

@test norm(sol .- sol2) ≈ 0 atol=1e-10

sol = solve(prob,  Rodas4())
sol2 = solve(prob2,Rodas4())

@test norm(sol .- sol2) ≈ 0 atol=1e-11

sol = solve(prob,  Rodas5())
sol2 = solve(prob2,Rodas5())

@test norm(sol .- sol2) ≈ 0 atol=1e-9

sol = solve(prob,  ImplicitEuler())
sol2 = solve(prob2,ImplicitEuler())

@test norm(sol .- sol2) ≈ 0 atol=1e-9

#sol = solve(prob, SDIRK2())
#sol2 = solve(prob2, SDIRK2())

#@test norm(sol .- sol2) ≈ 0 atol=1e-11



#sol = solve(prob,   Rodas3(),adaptive=false,dt=1/10)
#sol2 = solve(prob2, Rodas3(),adaptive=false,dt=1/10)
