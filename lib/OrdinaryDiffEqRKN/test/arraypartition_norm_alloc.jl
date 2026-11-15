using DiffEqBase, OrdinaryDiffEqRKN, Test

step_allocated_bytes(integrator) = @allocated step!(integrator)
residual_allocated_bytes(out, err, u₀, u₁) = @allocated DiffEqBase.calculate_residuals!(
    out, err, u₀, u₁, 1.0, 1.0, DiffEqBase.ODE_DEFAULT_NORM, 0.0
)

function kepler_acc_scalar!(ddu, du, u, p, t)
    r2 = u[1] * u[1] + u[2] * u[2]
    r3 = r2 * sqrt(r2)
    ddu[1] = -u[1] / r3
    ddu[2] = -u[2] / r3
    return nothing
end

@testset "ArrayPartition error norm allocations" begin
    prob = SecondOrderODEProblem(
        kepler_acc_scalar!, [0.0, 1.0], [1.0, 0.0], (0.0, 2π)
    )
    integrator = init(
        prob, DPRKN6(); reltol = 1.0e-8, abstol = 1.0e-8,
        save_everystep = false, dense = false
    )
    @test integrator.u isa DiffEqBase.RecursiveArrayTools.ArrayPartition
    step!(integrator)
    step!(integrator)
    step!(integrator)
    @test step_allocated_bytes(integrator) == 0
end

@testset "Heterogeneous ArrayPartition residual allocations" begin
    RAT = DiffEqBase.RecursiveArrayTools
    out = RAT.ArrayPartition(zeros(2), zeros(2, 2))
    err = RAT.ArrayPartition(ones(2), ones(2, 2))
    u₀ = RAT.ArrayPartition(ones(2), ones(2, 2))
    u₁ = RAT.ArrayPartition(fill(2.0, 2), fill(2.0, 2, 2))
    residual_allocated_bytes(out, err, u₀, u₁)
    @test residual_allocated_bytes(out, err, u₀, u₁) == 0
end
