using OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm,
    OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqImplicitAlgorithm,
    isadaptive, isdtchangeable, isimplicit
import OrdinaryDiffEqCore: alg_order, isfsal
using Test

struct GenericExplicitAlgorithm <: OrdinaryDiffEqAlgorithm end
struct GenericAdaptiveAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
struct GenericImplicitAlgorithm <: OrdinaryDiffEqImplicitAlgorithm end

alg_order(::GenericExplicitAlgorithm) = 1
alg_order(::GenericAdaptiveAlgorithm) = 4
alg_order(::GenericImplicitAlgorithm) = 2
isfsal(::GenericExplicitAlgorithm) = false

@testset "Generic algorithm trait contract" begin
    explicit = GenericExplicitAlgorithm()
    adaptive = GenericAdaptiveAlgorithm()
    implicit = GenericImplicitAlgorithm()

    @test alg_order(explicit) == 1
    @test alg_order(adaptive) == 4
    @test alg_order(implicit) == 2
    @test !isadaptive(explicit)
    @test isadaptive(adaptive)
    @test !isimplicit(explicit)
    @test isimplicit(implicit)
    @test !isfsal(explicit)
    @test isfsal(adaptive)
    @test isdtchangeable(implicit)
end
