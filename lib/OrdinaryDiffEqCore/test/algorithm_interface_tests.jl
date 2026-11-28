using OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm,
    OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqImplicitAlgorithm,
    isadaptive, isdtchangeable, isimplicit, CompositeAlgorithm,
    get_current_resumes_dt_after_tstop
import OrdinaryDiffEqCore: alg_order, isfsal, resumes_dt_after_tstop
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

struct ResumingAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
struct NonResumingAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
resumes_dt_after_tstop(::NonResumingAlgorithm) = false

@testset "tstop resumption follows the active composite member" begin
    @test get_current_resumes_dt_after_tstop(ResumingAlgorithm(), nothing)
    @test !get_current_resumes_dt_after_tstop(NonResumingAlgorithm(), nothing)
    for (algs, first_resumes) in (
            ((ResumingAlgorithm(), NonResumingAlgorithm()), true),
            ((NonResumingAlgorithm(), ResumingAlgorithm()), false),
        )
        composite = CompositeAlgorithm(algs, integrator -> 1)
        @test get_current_resumes_dt_after_tstop(composite, (; current = 1)) == first_resumes
        @test get_current_resumes_dt_after_tstop(composite, (; current = 2)) == !first_resumes
    end
end
