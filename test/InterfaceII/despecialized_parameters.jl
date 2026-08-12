using OrdinaryDiffEq, Test
import DiffEqBase
import SciMLBase

@test OrdinaryDiffEq.AutoDespecialize === SciMLBase.AutoDespecialize
@test OrdinaryDiffEq.AutoRespecialize === SciMLBase.AutoRespecialize
@test OrdinaryDiffEq.AutoDePSpecialize === SciMLBase.AutoDePSpecialize

struct FirstParameterLayout
    rate::Float64
end

struct SecondParameterLayout
    rate::Float64
    unused::Int
end

function despecialized_decay!(du, u, p, t)
    return du[1] = -p.rate * u[1]
end

function despecialized_problem(p)
    return ODEProblem{true, SciMLBase.AutoDespecialize}(
        despecialized_decay!, [1.0], (0.0, 1.0), p
    )
end

first_prob = despecialized_problem(FirstParameterLayout(1.0))
second_prob = despecialized_problem(SecondParameterLayout(2.0, 1))

@test first_prob.p isa FirstParameterLayout
@test second_prob.p isa SecondParameterLayout

first_concrete = DiffEqBase.get_concrete_problem(first_prob, true; alg = Tsit5())
second_concrete = DiffEqBase.get_concrete_problem(second_prob, true; alg = Tsit5())
@test typeof(first_concrete.p) === typeof(second_concrete.p) ===
    SciMLBase.DespecializedParameters
@test typeof(first_concrete) === typeof(second_concrete)
@test solve(
    first_prob, Tsit5(); saveat = [1.0], abstol = 1.0e-10, reltol = 1.0e-10
).u[end][1] ≈ exp(-1)
@test solve(
    second_prob, Tsit5(); saveat = [1.0], abstol = 1.0e-10, reltol = 1.0e-10
).u[end][1] ≈ exp(-2)
@test solve(
    first_prob, Rosenbrock23(); saveat = [1.0], abstol = 1.0e-12, reltol = 1.0e-12
).u[end][1] ≈ exp(-1)
