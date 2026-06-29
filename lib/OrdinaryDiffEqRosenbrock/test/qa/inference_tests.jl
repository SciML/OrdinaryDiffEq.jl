using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
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
# second dispatch. Both are now fixed, so `perform_step!` is free of runtime
# dispatch within OrdinaryDiffEqRosenbrock.
@testset "Rosenbrock perform_step! Inference Tests" begin
    function simple_system!(du, u, p, t)
        du .= t .* u
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise, and a
    # 50-dim system so the default linear solver takes its generic dispatch path.
    prob = ODEProblem{true, FullSpecialize}(simple_system!, ones(50), (0.0, 1.0))

    rosenbrock_solvers = [
        Rosenbrock23(), Rosenbrock32(), RosShamp4(), Veldd4(), Velds4(), GRK4T(), GRK4A(),
        Rodas3(), Rodas23W(), Rodas3P(), Rodas4(), Rodas42(), Rodas4P(), Rodas4P2(), Rodas5(),
        Rodas5P(), Rodas5Pe(), Rodas5Pr(), Rodas6P(),
    ]

    for solver in rosenbrock_solvers
        @testset "$(typeof(solver).name.name) perform_step! inference" begin
            integrator = init(
                prob, solver, dt = 0.1, save_everystep = false,
                abstol = 1.0e-10, reltol = 1.0e-10
            )
            step!(integrator)
            cache = integrator.cache

            JET.@test_opt target_modules = (OrdinaryDiffEqRosenbrock,) OrdinaryDiffEqCore.perform_step!(
                integrator, cache
            )
        end
    end
end
