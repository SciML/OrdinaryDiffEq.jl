using SciMLTesting, OrdinaryDiffEqVerner, Test

run_qa(
    OrdinaryDiffEqVerner;
    explicit_imports = true,
    ei_kwargs = (
        # `SciMLBase` module name is brought into scope by `@reexport using SciMLBase`
        # (intentional: this package re-exports the SciMLBase public API and also uses
        # it as a qualified accessor, e.g. `SciMLBase.solve`, `SciMLBase.AutoSpecialize`).
        no_implicit_imports = (; ignore = (:SciMLBase,)),
        # Qualified accesses to names that are not (yet) declared public by their owners.
        all_qualified_accesses_are_public = (;
            ignore = (
                # SciMLBase internals
                :FunctionWrapperSpecialize, :NoSpecialize, :has_lazy_interpolation,
                # OrdinaryDiffEqCore internals
                :increment_nf!, :set_EEst!, :lorenz, :lorenz_oop,
            )),
        # Explicit imports of names not (yet) declared public by their owners.
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internals (solver/cache/interpolation interface)
                :accept_step_controller, :alg_cache, :alg_stability_size, :AutoAlgSwitch,
                Symbol("@cache"), :CompiledFloats, :CompositeAlgorithm, :constvalue,
                :DerivativeOrderNotPossibleError, :explicit_rk_docstring, Symbol("@fold"),
                :get_fsalfirstlast, :isfsal, :_ode_addsteps!, :_ode_interpolant,
                :_ode_interpolant!, Symbol("@OnDemandTableauExtract"),
                :OrdinaryDiffEqAdaptiveAlgorithm, :OrdinaryDiffEqConstantCache,
                :OrdinaryDiffEqMutableCache, :perform_step!, :trivial_limiter!,
                :unwrap_alg,
                # SciMLBase internals
                Symbol("@def"), :_unwrap_val,
                # DiffEqBase internals
                :calculate_residuals, :calculate_residuals!,
                # TruncatedStacktraces internal macro
                Symbol("@truncate_stacktrace"),
            )),
    ),
)
