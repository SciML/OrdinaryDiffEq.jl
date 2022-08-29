using OrdinaryDiffEq, Test
function f(du, u, p, t)
    du[1] = 0.2u[1]
    du[2] = 0.4u[2]
end
u0 = ones(2)
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan, Float64[])

function lorenz(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
lorenzprob = ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[])
@test typeof(prob) === typeof(lorenzprob)
@test prob.f.f isa SciMLBase.FunctionWrappersWrappers.FunctionWrappersWrapper

t1 = @elapsed sol = solve(lorenzprob, Rosenbrock23())
t2 = @elapsed sol = solve(lorenzprob, Rosenbrock23(autodiff = false))

lorenzprob2 = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[])

t3 = @elapsed sol = solve(lorenzprob2, Rosenbrock23())
t4 = @elapsed sol = solve(lorenzprob2, Rosenbrock23(autodiff = false))

if VERSION >= v"1.8"
    @test 3t1 < t3
    @test t2 < t4
end

solve(prob, EPIRK4s3A(), dt = 1e-1)

function f_oop(u, p, t)
    [0.2u[1], 0.4u[2]]
end
u0 = ones(2)
tspan = (0.0, 1.0)
prob = ODEProblem{false, SciMLBase.FunctionWrapperSpecialize}(f_oop, u0, tspan, Float64[])

function lorenz_oop(u, p, t)
    [10.0(u[2] - u[1]), u[1] * (28.0 - u[3]) - u[2], u[1] * u[2] - (8 / 3) * u[3]]
end
lorenzprob = ODEProblem{false, SciMLBase.FunctionWrapperSpecialize}(lorenz_oop, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[])
@test typeof(prob) === typeof(lorenzprob)

# This one is fundamentally hard / broken
# Since the equation is not dependent on `t`, the output is not dual of t
# This is problem-dependent, so it is hard to deduce a priori
@test_broken t1 = @elapsed sol = solve(lorenzprob, Rosenbrock23())

t2 = @elapsed sol = solve(lorenzprob, Rosenbrock23(autodiff = false))

lorenzprob2 = ODEProblem{false, SciMLBase.NoSpecialize}(lorenz_oop, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[])

t3 = @elapsed sol = solve(lorenzprob2, Rosenbrock23())
t4 = @elapsed sol = solve(lorenzprob2, Rosenbrock23(autodiff = false))

#@test 5t1 < t3
#@test t2 < t4
