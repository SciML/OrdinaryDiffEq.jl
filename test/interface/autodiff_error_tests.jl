using OrdinaryDiffEq, Test, ADTypes
using OrdinaryDiffEqDifferentiation

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
@test_throws OrdinaryDiffEqDifferentiation.FirstAutodiffJacError solve(prob, Rosenbrock23())

function lorenz(u, p, t)
    du1 = 10.0(u[2] - u[1])
    a[1] = t
    du2 = u[1] * (28.0 - u[3]) - u[2]
    du3 = u[1] * u[2] - (8 / 3) * u[3]
    [du1, du2, du3]
end
@test_throws OrdinaryDiffEqDifferentiation.FirstAutodiffTgradError solve(
    prob, Rosenbrock23())

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    a[1] = u[2]
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 1.0)
prob = ODEProblem(lorenz!, u0, tspan)
@test_throws OrdinaryDiffEqDifferentiation.FirstAutodiffJacError solve(prob, Rosenbrock23())

function lorenz2!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    a[1] = t
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
prob = ODEProblem(lorenz2!, u0, tspan)
@test_throws OrdinaryDiffEqDifferentiation.FirstAutodiffTgradError solve(
    prob, Rosenbrock23())

## Test that nothing is using duals when autodiff=false
## https://discourse.julialang.org/t/rodas4-using-dual-number-for-time-with-autodiff-false/98256

for alg in [
    Rosenbrock23(autodiff = AutoFiniteDiff()),
    Rodas4(autodiff = AutoFiniteDiff()),
    Rodas5(autodiff = AutoFiniteDiff()),
    QNDF(autodiff = AutoFiniteDiff()),
    TRBDF2(autodiff = AutoFiniteDiff()),
    KenCarp4(autodiff = AutoFiniteDiff())
]
    u = [0.0, 0.0]
    function f1(u, p, t)
        #display(typeof(t))
        du = zeros(2)
        du[1] = 0.1 * u[1] + 0.2 * u[2]
        du[2] = 0.1 * t
        return du
    end
    prob = ODEProblem(f1, u, (0.0, 1.0))
    sol = solve(prob, alg)

    function f2(du, u, p, t)
        #display(typeof(t))
        du2 = zeros(2)
        du2[1] = 0.1 * u[1] + 0.2 * u[2]
        du2[2] = 0.1 * t
        du .= du2
    end
    prob = ODEProblem(f2, u, (0.0, 1.0))
    sol = solve(prob, alg)
end
