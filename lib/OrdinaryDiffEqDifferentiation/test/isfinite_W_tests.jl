using OrdinaryDiffEqDifferentiation: isfinite_W
using Test

struct OnlyIsfinite <: Number
    x::Float64
end
Base.isfinite(a::OnlyIsfinite) = isfinite(a.x)

# A dense matrix without scalar indexing, standing in for GPU arrays (which are
# `StridedMatrix` too): whole-array reductions work, elementwise iteration throws.
struct NoScalarMatrix <: DenseMatrix{Float64}
    data::Matrix{Float64}
end
Base.size(A::NoScalarMatrix) = size(A.data)
Base.getindex(::NoScalarMatrix, ::Int...) = error("scalar indexing")
Base.all(f::typeof(isfinite), A::NoScalarMatrix) = all(f, A.data)

@testset "isfinite_W" begin
    for T in (Float16, Float32, Float64, ComplexF32, ComplexF64, BigFloat)
        W = T[1 2; 3 4]
        @test isfinite_W(W)
        @test isfinite_W(view(W, :, 1:2))
        for bad in (Inf, -Inf, NaN)
            Wbad = copy(W)
            Wbad[2, 1] = bad
            @test !isfinite_W(Wbad)
            @test !isfinite_W(view(Wbad, 1:2, 1:2))
        end
    end
    W = ComplexF64[1 2; 3 4]
    W[1, 2] = complex(1.0, Inf)
    @test !isfinite_W(W)

    @test isfinite_W([1 // 2 1; 0 3])
    @test !isfinite_W([1 // 0 1; 0 3])
    @test !isfinite_W(fill(-1 // 0, 2, 2))
    @test isfinite_W([1 2; 3 4])

    @test isfinite_W(fill(OnlyIsfinite(1.0), 2, 2))
    @test !isfinite_W([OnlyIsfinite(1.0) OnlyIsfinite(NaN)])

    @test isfinite_W(NoScalarMatrix([1.0 2.0; 3.0 4.0]))
    @test !isfinite_W(NoScalarMatrix([1.0 NaN; 3.0 4.0]))

    @test isfinite_W(2.0)
    @test !isfinite_W(Inf)
    @test isfinite_W(nothing)
end
