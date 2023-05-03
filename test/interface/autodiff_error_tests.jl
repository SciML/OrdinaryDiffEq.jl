using OrdinaryDiffEq, Test

const a = Float64[1.0]

function lorenz(u, p, t)
    du1 = 10.0(u[2] - u[1])
    a[1] = u[2]
    du2 = u[1] * (28.0 - u[3]) - u[2]
    du3 = u[1] * u[2] - (8 / 3) * u[3]
    [du1, du2, du3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 1.0)
prob = ODEProblem(lorenz, u0, tspan)
@test_throws OrdinaryDiffEq.FirstAutodiffJacError solve(prob, Rosenbrock23())

function lorenz(u, p, t)
    du1 = 10.0(u[2] - u[1])
    a[1] = t
    du2 = u[1] * (28.0 - u[3]) - u[2]
    du3 = u[1] * u[2] - (8 / 3) * u[3]
    [du1, du2, du3]
end
@test_throws OrdinaryDiffEq.FirstAutodiffTgradError solve(prob, Rosenbrock23())

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    a[1] = u[2]
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 1.0)
prob = ODEProblem(lorenz!, u0, tspan)
@test_throws OrdinaryDiffEq.FirstAutodiffJacError solve(prob, Rosenbrock23())

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    a[1] = t
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
@test_throws OrdinaryDiffEq.FirstAutodiffTgradError solve(prob, Rosenbrock23())

## Test that nothing is using duals when autodiff=false
## https://discourse.julialang.org/t/rodas4-using-dual-number-for-time-with-autodiff-false/98256

for alg in [Rosenbrock23(autodiff=false), Rodas4(autodiff=false), Rodas5(autodiff=false), QNDF(autodiff=false), TRBDF2(autodiff=false), KenCarp4(autodiff=false)]
    u = [0.0, 0.0]
    function du(u,p,t)
        #display(typeof(t))
        du = zeros(2)
        du[1] = 0.1*u[1] + 0.2*u[2]
        du[2] = 0.1*t
        return du
    end
    prob = ODEProblem(du,u,(0.0,1.0))
    sol = solve(prob,alg)
end