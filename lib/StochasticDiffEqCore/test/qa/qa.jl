using SciMLTesting, StochasticDiffEqCore, Test
using JET

# StochasticDiffEqCore is an extension package for OrdinaryDiffEqCore's solver
# loop: it dispatches on and extends OrdinaryDiffEqCore-internal traits/loop
# helpers and a handful of still-internal SciMLBase/DiffEqBase interface names.
# Every remaining EI ignore below is a genuinely non-public (or non-owner) name
# from an upstream package; each is documented with its source. Names that
# became public (SciMLBase abstract types, is_diagonal_noise, __solve/__init,
# DEVerbosity, etc.) were migrated to their public owners in src/ rather than
# ignored. See SciML/OrdinaryDiffEq.jl#3776.

# OrdinaryDiffEqCore-internal names that this package dispatches on or extends.
# These are not part of OrdinaryDiffEqCore's public API; they need make-public
# upstream before they can be dropped from the ignore lists.
const ODEC_INTERNAL = (
    # Algorithm/cache supertypes + controllers/defaults (explicit imports)
    :ODEIntegrator, :StochasticDiffEqAlgorithm, :StochasticDiffEqAdaptiveAlgorithm,
    :StochasticDiffEqCompositeAlgorithm, :StochasticDiffEqRODEAlgorithm,
    :StochasticDiffEqRODEAdaptiveAlgorithm, :StochasticDiffEqRODECompositeAlgorithm,
    :StochasticDiffEqNewtonAdaptiveAlgorithm, :StochasticDiffEqNewtonAlgorithm,
    :StochasticDiffEqJumpAlgorithm, :StochasticDiffEqJumpAdaptiveAlgorithm,
    :StochasticDiffEqJumpNewtonAdaptiveAlgorithm, :StochasticDiffEqJumpDiffusionAlgorithm,
    :StochasticDiffEqJumpDiffusionAdaptiveAlgorithm,
    :StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm,
    :StochasticDiffEqCache, :StochasticDiffEqConstantCache, :StochasticDiffEqMutableCache,
    :beta1_default, :beta2_default,
    :qmin_default, :qmax_default, :qsteady_min_default, :qsteady_max_default,
    :issplit, :is_composite_algorithm, :perform_step!, :handle_callback_modifiers!,
    # Loop/trait/autodiff helpers extended via qualified `OrdinaryDiffEqCore.x`
    :DEOptions, :_determine_initdt, :_get_fdtype, :_get_fwd_chunksize,
    :_get_fwd_chunksize_int, :_initialize_dae!, :_ode_init, :accept_noise!,
    :alg_autodiff, :alg_difftype, :alg_extrapolates, :concrete_jac, :gamma_default,
    :get_chunksize, :get_current_alg_autodiff, :get_fsalfirstlast, :has_autodiff,
    :is_composite_cache, :is_constant_cache, :is_noise_saveable, :isfsal,
    :noise_curt, :ode_determine_initdt, :reinit_noise!, :reject_noise!,
    :save_noise!, :standardtag,
)

# Still-internal SciMLBase interface names (problem/alg supertypes, alias
# specifiers, mass-matrix/initialization helpers, autodiff val unwrap).
const SCIMLBASE_INTERNAL = (
    :AbstractRODEProblem, :AbstractSDDEProblem, :AbstractSDDEIntegrator,
    :AlgorithmInterpretation, :alg_interpretation, :RODEAliasSpecifier,
    :SDEAliasSpecifier, :__has_mass_matrix, :has_initializeprob,
    :parameterless_type, :_unwrap_val,
)

# Still-internal DiffEqBase names. `@..` is owned by FastBroadcast and obtained
# through DiffEqBase's re-export (the standard SciML access path), so it is also
# listed in the via-owners ignore below.
const DIFFEQBASE_INTERNAL = (:prob2dtmin, :ODE_DEFAULT_UNSTABLE_CHECK, Symbol("@.."))

# Non-public names from other upstream packages.
const JUMPPROCESSES_INTERNAL = (:reset_jump_problem!, :resetted_jump_problem)
const DIFFEQNOISEPROCESS_INTERNAL = (:resize_stack!,)
const FORWARDDIFF_INTERNAL = (:Tag, :pickchunksize)
const BASE_INTERNAL = (Symbol("@pure"),)

const NONPUBLIC_IGNORE = (
    ODEC_INTERNAL..., SCIMLBASE_INTERNAL..., DIFFEQBASE_INTERNAL...,
    JUMPPROCESSES_INTERNAL..., DIFFEQNOISEPROCESS_INTERNAL...,
    FORWARDDIFF_INTERNAL..., BASE_INTERNAL...,
)

run_qa(
    StochasticDiffEqCore;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
    ei_kwargs = (;
        # `@..` reaches StochasticDiffEqCore via DiffEqBase's re-export of
        # FastBroadcast; FastBroadcast is not a direct dependency.
        all_explicit_imports_via_owners = (; ignore = (Symbol("@.."),)),
        all_qualified_accesses_are_public = (; ignore = NONPUBLIC_IGNORE),
        all_explicit_imports_are_public = (; ignore = NONPUBLIC_IGNORE),
    ),
)
