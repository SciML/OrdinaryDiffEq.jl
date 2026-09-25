using OrdinaryDiffEqCore
using Test

# OrdinaryDiffEq.jl#1402 / JuliaPy/PyCall.jl#900:
# some wrappers return Vector{Bool} from `==`, which is invalid in boolean context.
struct ArrayValuedEq
    x::Vector{Float64}
end
Base.:(==)(a::ArrayValuedEq, b::ArrayValuedEq) = a.x .== b.x

@testset "_bool_equal coerces array-valued ==" begin
    @test OrdinaryDiffEqCore._bool_equal(1.0, 1.0)
    @test !OrdinaryDiffEqCore._bool_equal(1.0, 2.0)
    @test OrdinaryDiffEqCore._bool_equal([1.0, 2.0], [1.0, 2.0])
    @test !OrdinaryDiffEqCore._bool_equal([1.0, 2.0], [1.0, 3.0])

    a = ArrayValuedEq([1.0, 2.0, 3.0])
    b = ArrayValuedEq([1.0, 2.0, 3.0])
    c = ArrayValuedEq([1.0, 2.0, 4.0])
    @test (a == b) isa AbstractArray
    @test OrdinaryDiffEqCore._bool_equal(a, b)
    @test !OrdinaryDiffEqCore._bool_equal(a, c)
    # boolean context must work
    @test OrdinaryDiffEqCore._bool_equal(a, b) && true
end
