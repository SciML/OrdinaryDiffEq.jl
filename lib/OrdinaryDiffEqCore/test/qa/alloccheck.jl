using OrdinaryDiffEqCore
using Test

# Core infrastructure function tests (loopheader!, _loopfooter!, apply_step!, etc.)
# require an ODEIntegrator, which requires a concrete solver subpackage.
# All solver subpackages depend on OrdinaryDiffEqCore, creating a circular test
# dependency that Julia's Pkg sandbox cannot resolve.
# These tests live in each solver's own allocation_tests.jl instead.

@testset "OrdinaryDiffEqCore AllocCheck" begin
    @test true
end
