using Test, RecursiveArrayTools
using DiffEqBase: calculate_residuals!, ODE_DEFAULT_NORM

expand_tolerance(tol::Number, n) = fill(tol, n)
expand_tolerance(tol, n) = vec(collect(tol))

function reference_residuals(ũ, u₀, u₁, α, ρ, t)
    n = length(ũ)
    αs = expand_tolerance(α, n)
    ρs = expand_tolerance(ρ, n)
    return [
        ũ[i] / (αs[i] + max(ODE_DEFAULT_NORM(u₀[i], t), ODE_DEFAULT_NORM(u₁[i], t)) * ρs[i])
            for i in 1:n
    ]
end

@testset "ArrayPartition residuals ($layout): $(nameof(typeof(α))) abstol, $(nameof(typeof(ρ))) reltol" for (layout, partition) in (
            ("Vector, Matrix", v -> ArrayPartition(v[1:2], reshape(v[3:6], 2, 2))),
            ("Vector, Vector", v -> ArrayPartition(v[1:3], v[4:6])),
        ),
        α in (1.0e-3, partition(collect(1.0:6.0)), collect(1.0:6.0)),
        ρ in (0.5, partition(collect(0.1:0.1:0.6)), collect(0.1:0.1:0.6))

    ũ = partition([0.5, -1.0, 2.0, 1.5, -0.25, 3.0])
    u₀ = partition([1.0, 2.0, -3.0, 0.5, 4.0, -1.0])
    u₁ = partition([-2.0, 1.0, 1.0, 2.0, -5.0, 0.25])
    out = zero(ũ)
    calculate_residuals!(out, ũ, u₀, u₁, α, ρ, ODE_DEFAULT_NORM, 0.0)
    @test collect(out) ≈ reference_residuals(ũ, u₀, u₁, α, ρ, 0.0)
end
