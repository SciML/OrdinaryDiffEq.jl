using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using LinearAlgebra
using JET
using Test

# Regression test for https://github.com/SciML/OrdinaryDiffEq.jl/issues/3800.
#
# `perform_step!` passed the system matrix to the linear solve through a ternary
# keyword argument, `A = (repeat_step || !new_W) ? nothing : W`. That makes `A`
# inferred as `Union{Nothing, typeof(W)}`, which turns the `dolinsolve` keyword
# call into a runtime-dispatched call inside the otherwise type-stable
# `perform_step!`. A related dead branch (`reshape(cache.algebraic_vars, ...)`
# guarded on `mass_matrix !== I` instead of `algebraic_vars !== nothing`) added a
# second dispatch in the non-DAE case. Both are now fixed, so `perform_step!` is
# free of runtime dispatch within OrdinaryDiffEqRosenbrock.
#
# Both an ODE (`mass_matrix === I`, `algebraic_vars === nothing`) and a mass-matrix
# DAE (`algebraic_vars isa Vector`, which exercises the `reshape`/`ifelse` branch)
# are checked, since they take different paths through `perform_step!`.
@testset "Rosenbrock perform_step! Inference Tests" begin
    rosenbrock_solvers = [
        Rosenbrock23(), Rosenbrock32(), RosShamp4(), Veldd4(), Velds4(), GRK4T(), GRK4A(),
        Rodas3(), Rodas3d(), Rodas23W(), Rodas3P(), Rodas4(), Rodas42(), Rodas4P(), Rodas4P2(), Rodas5(),
        Rodas5P(), Rodas5Pe(), Rodas5Pr(), Rodas6P(),
    ]

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise, and a
    # 50-dim system so the default linear solver takes its generic dispatch path.
    function ode_system!(du, u, p, t)
        du .= t .* u
    end
    ode_prob = ODEProblem{true, FullSpecialize}(ode_system!, ones(50), (0.0, 1.0))

    # Index-1 mass-matrix DAE: u[i]' = -u[i] for the first half, and the algebraic
    # constraint u[n+i] = u[i] for the second half (M = diag(ones(n), zeros(n))).
    function dae_system!(du, u, p, t)
        n = length(u) ÷ 2
        for i in 1:n
            du[i] = -u[i]
            du[n + i] = u[i] - u[n + i]
        end
    end
    n = 25
    mass_matrix = Matrix(Diagonal([ones(n); zeros(n)]))
    dae_prob = ODEProblem(
        ODEFunction{true, FullSpecialize}(dae_system!; mass_matrix = mass_matrix),
        ones(2n), (0.0, 1.0)
    )

    for (label, prob, tol) in (
            ("ODE", ode_prob, 1.0e-10),
            ("DAE", dae_prob, 1.0e-8),
        )
        @testset "$label" begin
            for solver in rosenbrock_solvers
                @testset "$(typeof(solver).name.name) perform_step! inference" begin
                    integrator = init(
                        prob, solver, dt = 0.05, save_everystep = false,
                        abstol = tol, reltol = tol
                    )
                    step!(integrator)
                    cache = integrator.cache

                    JET.@test_opt target_modules = (OrdinaryDiffEqRosenbrock,) OrdinaryDiffEqCore.perform_step!(
                        integrator, cache
                    )
                end
            end
        end
    end
end
