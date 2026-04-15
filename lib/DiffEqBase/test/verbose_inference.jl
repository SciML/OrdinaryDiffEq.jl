using DiffEqBase, Test

# Verify that the Bool path of `_process_verbose_param` is type-stable when
# the Bool is a compile-time literal (the default `verbose = true` case from
# `solve`/`init`). Without `Base.@constprop :aggressive`, both branches return
# different concrete `DEVerbosity{...}` parameterizations and the result is a
# `Union{...}` that poisons inference all the way through `solve`.

@testset "_process_verbose_param inference" begin
    # Val-dispatch is unconditionally type-stable.
    @inferred DiffEqBase._process_verbose_param(Val(true))
    @inferred DiffEqBase._process_verbose_param(Val(false))

    # Constprop on literal Bool resolves to a single concrete type.
    f_true() = DiffEqBase._process_verbose_param(true)
    f_false() = DiffEqBase._process_verbose_param(false)
    @test Base.return_types(f_true, ())[1] === typeof(DiffEqBase.DEFAULT_VERBOSE)
    @test Base.return_types(f_false, ())[1] === typeof(DiffEqBase.NONE_VERBOSE)
end
