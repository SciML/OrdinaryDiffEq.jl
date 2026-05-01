using DiffEqBase, Test

# Verify that the DEVerbosity paths of `_process_verbose_param` are type-stable
# and that passing a Bool throws (v7 breaking change).

@testset "_process_verbose_param inference" begin
    # DEVerbosity passthrough is type-stable.
    @inferred DiffEqBase._process_verbose_param(DiffEqBase.DEFAULT_VERBOSE)

    # Preset dispatch is type-stable.
    @inferred DiffEqBase._process_verbose_param(DiffEqBase.SciMLLogging.Standard())

    # Bool is no longer accepted — throws ArgumentError.
    @test_throws ArgumentError DiffEqBase._process_verbose_param(true)
    @test_throws ArgumentError DiffEqBase._process_verbose_param(false)
end
