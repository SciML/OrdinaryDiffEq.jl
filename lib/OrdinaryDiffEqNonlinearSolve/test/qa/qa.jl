using SciMLTesting, OrdinaryDiffEqNonlinearSolve, Test
using Aqua

@static if VERSION >= v"1.11.0-DEV.469"
    for name in (
            :NLNewton, :NLFunctional, :NLAnderson, :NonlinearSolveAlg,
            :HomotopyNonlinearSolveAlg, :build_nlsolver, :nlsolve!, :nlsolvefail,
            :markfirststage!, :du_alias_or_new, :can_smooth_est, :compute_step!,
            :initial_η, :anderson, :anderson!,
        )
        @test Base.ispublic(OrdinaryDiffEqNonlinearSolve, name)
    end
    for name in (
            :NLSolver, :NLFunctionalCache, :NLAndersonCache, :NLNewtonCache,
            :NonlinearSolveCache, :HomotopyNonlinearSolveCache,
        )
        @test !Base.ispublic(OrdinaryDiffEqNonlinearSolve, name)
    end
end

# Ambiguities are checked below over this package's own methods instead of Aqua's default
# `[pkg, Base, Core]` module set; see the `Aqua.test_ambiguities` call.
run_qa(
    OrdinaryDiffEqNonlinearSolve;
    # `BrownFullBasicInit`/`ShampineCollocationInit` are the user-facing DAE
    # initialization algorithms. DiffEqBase owns, exports, and documents them; this
    # package implements their solver-side methods and re-exports the names so
    # `using OrdinaryDiffEq` keeps surfacing them unqualified.
    reexports_allow = (:BrownFullBasicInit, :ShampineCollocationInit),
    aqua_kwargs = (; piracies = false, ambiguities = false),
    ei_kwargs = (;
        # Imported into this namespace for dependent sublibraries; ExplicitImports
        # cannot see that cross-package use and reports them as stale.
        no_stale_explicit_imports = (; ignore = (:MatrixOperator, :_reshape)),
        # Names imported from a re-exporter rather than their defining package.
        # `@SciMLMessage` is owned by SciMLLogging but reached through OrdinaryDiffEqCore
        # (SciMLLogging is not a direct dependency); `WOperator`/`StaticWOperator` are the
        # OrdinaryDiffEqDifferentiation W-operator types, which ExplicitImports attributes
        # to SciMLOperators (their abstract supertype's package).
        all_explicit_imports_via_owners = (;
            ignore = (Symbol("@SciMLMessage"), :WOperator, :StaticWOperator),
        ),
        # Non-public names of upstream packages with no public alternative. Each is owned
        # by the package it is imported from; remove an entry once that package marks the
        # name `public`/exports it. Tracked in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                # Base — owner-internal, no public alternative
                :RefValue,
                # NonlinearSolveBase — owner-internal, no public alternative
                :not_terminated, :get_u, :get_fu,
                # OrdinaryDiffEqDifferentiation — owner-internal, no public alternative
                :default_krylov_warm_start,
                # `@SciMLMessage` reached through OrdinaryDiffEqCore (owner SciMLLogging).
                Symbol("@SciMLMessage"),
                # SciMLBase internals (owner of these names but not public).
                :_vec, :_reshape, :postamble!,
                # SciMLOperators abstract type (not public).
                :AbstractSciMLOperator,
                # OrdinaryDiffEqDifferentiation W-operator types (attributed to SciMLOperators).
                :WOperator, :StaticWOperator,
                # OrdinaryDiffEqDifferentiation:
                :build_J_W, :build_jac_config, :build_uf, :dolinsolve,
                :is_always_new, :jacobian!, :jacobian2W!, :resize_jac_config!,
                :update_W!, :wrapprecs,
                # ForwardDiff / StaticArraysCore internals.
                :Dual, :StaticArray,
            ),
        ),
        # Qualified `Module.name` accesses to non-public names with no public alternative.
        all_qualified_accesses_are_public = (;
            ignore = (
                # Base — owner-internal, no public alternative
                :RefValue,
                # NonlinearSolveBase — owner-internal, no public alternative
                :not_terminated,
                # OrdinaryDiffEqDifferentiation — owner-internal, no public alternative
                :default_krylov_warm_start,
                # SciMLBase internals.
                :value, :anyeltypedual, :forwarddiff_chunksize,
                :has_Wfact, :has_Wfact_t, :has_jac_u, :has_jac_du, :postamble!,
                # ForwardDiff internals.
                :Tag, :pickchunksize,
            ),
        ),
    ),
)

# `Aqua.test_all` scans `[pkg, Base, Core]`, which reports every ambiguity between a Base
# method and any loaded package — for this dependency tree, ForwardDiff's `Dual`
# constructors, `convert`, `==`, `^` and `log` against Base, none of which involve a method
# defined here or can be resolved from here. `detect_ambiguities` keeps only pairs with a
# method defined in the listed modules, so scanning this package alone still catches every
# ambiguity this package is responsible for — including ones against Base or a dependency —
# while dropping the ForwardDiff/Base noise.
@testset "Method ambiguity (own methods)" begin
    Aqua.test_ambiguities(OrdinaryDiffEqNonlinearSolve)
end
