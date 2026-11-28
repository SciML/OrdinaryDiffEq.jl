using OrdinaryDiffEqRosenbrock, ReverseDiff, Test
using ADTypes: AutoFiniteDiff
using DiffEqBase: SensitivityADPassThrough
using SciMLBase: successful_retcode

function coupled_rhs!(du, u, p, t)
    du[1] = -p[1] * u[1]
    du[2] = p[1] * u[1] - u[2]
    return nothing
end

function uncoupled_rhs!(du, u, p, t)
    du[1] = -2 * u[1]
    du[2] = -u[2]
    return nothing
end

function first_state(f, u0, p, tspan, alg)
    sol = solve(
        ODEProblem(f, u0, tspan, p), alg;
        sensealg = SensitivityADPassThrough(), abstol = 1.0e-8, reltol = 1.0e-8
    )
    @test successful_retcode(sol)
    return sol.u[end][1]
end

@testset "ReverseDiff through $(nameof(typeof(alg)))" for alg in (
        Rodas5(autodiff = AutoFiniteDiff()), Rosenbrock23(autodiff = AutoFiniteDiff()),
    )
    # u0 converted to eltype(p) gives W an eltype that `lu` does not preserve
    grad = ReverseDiff.gradient(
        p -> first_state(coupled_rhs!, eltype(p).([1.0, 0.0]), p, (0.0, 1.0), alg), [2.0]
    )
    @test grad ≈ [-exp(-2)] rtol = 1.0e-5

    grad = ReverseDiff.gradient(
        q -> first_state(uncoupled_rhs!, [1.0 + q[1], 1.0 + q[2]], nothing, (0.0, 1.0), alg),
        [1.0, 0.5]
    )
    @test grad ≈ [exp(-2), 0.0] rtol = 1.0e-5

    grad = ReverseDiff.gradient(
        q -> first_state(uncoupled_rhs!, [q[1], q[2]], nothing, (0.0, 0.0), alg), [1.0, 0.5]
    )
    @test grad == [1.0, 0.0]
end
