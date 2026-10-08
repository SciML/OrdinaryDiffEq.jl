using OrdinaryDiffEqDifferentiation: isfinite_W
using Test

struct OnlyIsfinite <: Number
    x::Float64
end
Base.isfinite(a::OnlyIsfinite) = isfinite(a.x)

@testset "isfinite_W" begin
    for T in (Float32, Float64, ComplexF64, BigFloat)
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

    @test isfinite_W(2.0)
    @test !isfinite_W(Inf)
    @test isfinite_W(nothing)
end
