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

@testset "global_error group" begin
    SL = DiffEqBase.SciMLLogging
    fields = (
        :global_error_event_mismatch, :global_error_companion_order, :global_error_final_solve,
        :global_error_ad_fallback, :global_error_grazing_event,
    )
    # The group keyword sets all of its fields, and a field keyword overrides the group.
    v = DEVerbosity(global_error = SL.Silent(), global_error_grazing_event = SL.WarnLevel())
    @test all(f -> getproperty(v, f) == SL.Silent(), fields[1:4])
    @test v.global_error_grazing_event == SL.WarnLevel()
    # None silences the fields; Standard warns for all of them.
    @test all(f -> getproperty(DEVerbosity(SL.None()), f) == SL.Silent(), fields)
    @test all(f -> getproperty(DEVerbosity(SL.Standard()), f) == SL.WarnLevel(), fields)
    @inferred DiffEqBase._process_verbose_param(v)
end
