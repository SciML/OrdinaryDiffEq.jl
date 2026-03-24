using DiffEqBase.InternalEuler, SciMLBase, DiffEqBase, Test

# Here's the problem to solve

struct LorenzFunction <: Function
    syms::Vector{Symbol}
end

function (::LorenzFunction)(u, p, t)
    return [10.0(u[2] - u[1]), u[1] * (28.0 - u[3]) - u[2], u[1] * u[2] - (8 / 3) * u[3]]
end
lorenz = LorenzFunction([:x, :y, :z])

u0 = [1.0, 5.0, 10.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)
dt = 0.1
sol = solve(prob, InternalEuler.FwdEulerAlg(), tstops = 0:dt:1)
syms = [:x, :y, :z]

@test SciMLBase.interpret_vars([(0, 1), (1, 3), (4, 5)], sol) == [
    (SciMLBase.DEFAULT_PLOT_FUNC, 0, 1),
    (SciMLBase.DEFAULT_PLOT_FUNC, 1, 3),
    (SciMLBase.DEFAULT_PLOT_FUNC, 4, 5),
]
@test SciMLBase.interpret_vars([1, (1, 3), (4, 5)], sol) == [
    (SciMLBase.DEFAULT_PLOT_FUNC, 0, 1),
    (SciMLBase.DEFAULT_PLOT_FUNC, 1, 3),
    (SciMLBase.DEFAULT_PLOT_FUNC, 4, 5),
]
@test SciMLBase.interpret_vars([1, 3, 4], sol) == [
    (SciMLBase.DEFAULT_PLOT_FUNC, 0, 1),
    (SciMLBase.DEFAULT_PLOT_FUNC, 0, 3),
    (SciMLBase.DEFAULT_PLOT_FUNC, 0, 4),
]
@test SciMLBase.interpret_vars(([1, 2, 3], [4, 5, 6]), sol) == [
    (SciMLBase.DEFAULT_PLOT_FUNC, 1, 4),
    (SciMLBase.DEFAULT_PLOT_FUNC, 2, 5),
    (SciMLBase.DEFAULT_PLOT_FUNC, 3, 6),
]
@test SciMLBase.interpret_vars((1, [2, 3, 4]), sol) == [
    (SciMLBase.DEFAULT_PLOT_FUNC, 1, 2),
    (SciMLBase.DEFAULT_PLOT_FUNC, 1, 3),
    (SciMLBase.DEFAULT_PLOT_FUNC, 1, 4),
]

f(x, y) = (x + y, y)
@test SciMLBase.interpret_vars([(f, 0, 1), (1, 3), (4, 5)], sol) ==
    [(f, 0, 1), (SciMLBase.DEFAULT_PLOT_FUNC, 1, 3), (SciMLBase.DEFAULT_PLOT_FUNC, 4, 5)]
@test SciMLBase.interpret_vars([1, (f, 1, 3), (4, 5)], sol) ==
    [(SciMLBase.DEFAULT_PLOT_FUNC, 0, 1), (f, 1, 3), (SciMLBase.DEFAULT_PLOT_FUNC, 4, 5)]
@test SciMLBase.interpret_vars([1, (f, 0, 1), (1, 2)], sol) ==
    [(SciMLBase.DEFAULT_PLOT_FUNC, 0, 1), (f, 0, 1), (SciMLBase.DEFAULT_PLOT_FUNC, 1, 2)]
@test SciMLBase.interpret_vars([(1, 2)], sol) ==
    [(SciMLBase.DEFAULT_PLOT_FUNC, 1, 2)]
@test SciMLBase.interpret_vars((f, 1, 2), sol) == [(f, 1, 2)]
