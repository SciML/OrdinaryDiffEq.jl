using StochasticDiffEqROCK
using Test

@testset "StochasticDiffEqROCK" begin
    @test SROCK1 <: StochasticDiffEqROCK.StochasticDiffEqAlgorithm
    @test SROCK2 <: StochasticDiffEqROCK.StochasticDiffEqAlgorithm
    @test KomBurSROCK2 <: StochasticDiffEqROCK.StochasticDiffEqAlgorithm
    @test SROCKC2 <: StochasticDiffEqROCK.StochasticDiffEqAlgorithm
    @test SROCKEM <: StochasticDiffEqROCK.StochasticDiffEqAlgorithm
    @test SKSROCK <: StochasticDiffEqROCK.StochasticDiffEqAlgorithm
    @test TangXiaoSROCK2 <: StochasticDiffEqROCK.StochasticDiffEqAlgorithm

    # Test constructors
    @test SROCK1() isa SROCK1
    @test SROCK2() isa SROCK2
    @test KomBurSROCK2() isa KomBurSROCK2
    @test SROCKC2() isa SROCKC2
    @test SROCKEM() isa SROCKEM
    @test SKSROCK() isa SKSROCK
    @test TangXiaoSROCK2() isa TangXiaoSROCK2
end
