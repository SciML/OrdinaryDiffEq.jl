using OrdinaryDiffEqCore
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqCore.

Note: Core infrastructure functions (loopheader!, _loopfooter!, apply_step!,
fix_dt_at_bounds!, modify_dt_for_tstops!, perform_step!) all require a concrete
ODEIntegrator, which requires importing a solver subpackage. Since all solver
subpackages depend on OrdinaryDiffEqCore, importing them here creates a circular
test dependency that Julia's Pkg sandbox cannot resolve.

These infrastructure function tests live in each solver's own allocation_tests.jl,
where the solver subpackage is the one under test (not Core).
See: OrdinaryDiffEqTsit5/test/allocation_tests.jl for loopheader!, _loopfooter!,
apply_step!, fix_dt_at_bounds!, modify_dt_for_tstops!, and perform_step! (Tsit5).

This file tests any Core-level functionality that can be verified without
importing a solver subpackage.
"""

@testset "OrdinaryDiffEqCore AllocCheck" begin
    # placeholder: infrastructure function tests are in solver-specific alloccheck files
    # (see OrdinaryDiffEqTsit5/test/allocation_tests.jl)
    @test true
end
