using OrdinaryDiffEqTaylorSeries, SparseArrays, Test

struct ReverseStorageVector{T} <: AbstractVector{T}
    storage::Vector{T}
end

Base.parent(u::ReverseStorageVector) = u.storage
Base.size(u::ReverseStorageVector) = size(parent(u))
Base.IndexStyle(::Type{<:ReverseStorageVector}) = IndexLinear()
Base.getindex(u::ReverseStorageVector, i::Int) = parent(u)[end + 1 - i]
Base.setindex!(u::ReverseStorageVector, x, i::Int) = (parent(u)[end + 1 - i] = x)
Base.copy(u::ReverseStorageVector) = ReverseStorageVector(copy(parent(u)))
Base.zero(u::ReverseStorageVector) = ReverseStorageVector(zero(parent(u)))
Base.convert(::Type{ReverseStorageVector{T}}, u::Vector{T}) where {T} = ReverseStorageVector(reverse(u))

function indexed_decay!(du, u, p, t)
    for i in eachindex(u, du)
        du[i] = -i * u[i]
    end
    return nothing
end

sparse_relaxation!(du, u, p, t) = (du .= p .- u; nothing)

@testset "In-place Taylor jet state indexing" begin
    @testset "$(nameof(typeof(solver)))" for solver in (
            ExplicitTaylor(order = Val(5)),
            ExplicitTaylorAdaptiveOrder(min_order = Val(3), max_order = Val(5)),
        )
        @testset "Logical order with reversed parent storage" begin
            u0 = ReverseStorageVector([1.0, 2.0, 3.0])
            prob = ODEProblem{true, FullSpecialize}(indexed_decay!, u0, (0.0, 0.01))
            sol = solve(prob, solver; dt = 0.01, adaptive = false)
            @test successful_retcode(sol)
            @test sol.u[end] isa ReverseStorageVector
            @test sol.u[end] ≈ [3.0, 2.0, 1.0] .* exp.(-0.01 .* (1:3)) rtol = 1.0e-9
        end

        @testset "Sparse initial state with implicit zeros" begin
            u0 = sparsevec([3], [1.0], 3)
            target = [1.0, 2.0, 3.0]
            prob = ODEProblem{true, FullSpecialize}(sparse_relaxation!, u0, (0.0, 0.1), target)
            sol = solve(prob, solver; dt = 0.005, adaptive = false)
            @test successful_retcode(sol)
            @test sol.u[end] isa SparseVector
            @test sol.u[end] ≈ target + (u0 - target) * exp(-0.1) rtol = 1.0e-9
        end
    end
end
