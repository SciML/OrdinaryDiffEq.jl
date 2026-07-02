using SciMLTesting, StochasticDiffEqCore, Test
using JET

# StochasticDiffEqCore extends OrdinaryDiffEqCore's solver loop and dispatches on
# a handful of still-internal SciMLBase/DiffEqBase interface names. The base
# packages (OrdinaryDiffEqCore/DiffEqBase) now declare their solver-author API
# public, so every name that became public was dropped from the ignore lists
# below. Each remaining ignore is a genuinely non-public (or non-owner) name
# from an upstream package, documented with its source. See
# SciML/OrdinaryDiffEq.jl#3776.

# OrdinaryDiffEqCore names this package dispatches on/extends that were kept
# owner-internal (loop/trait/autodiff/noise helpers not in the public
# solver-author surface). These need make-public upstream before they can drop.
const ODEC_INTERNAL = (
    :_determine_initdt, :_get_fdtype, :_get_fwd_chunksize, :_get_fwd_chunksize_int,
    :_initialize_dae!, :_ode_init, :accept_noise!, :concrete_jac, :get_chunksize,
    :handle_callback_modifiers!, :has_autodiff, :is_noise_saveable, :noise_curt,
    :ode_determine_initdt, :qsteady_max_default, :qsteady_min_default,
    :reinit_noise!, :reject_noise!, :save_noise!, :standardtag,
)

# Still-internal SciMLBase interface names (pending SciMLBase#1412 round-5
# make-public; not yet public on the registered SciMLBase this branch resolves).
const SCIMLBASE_INTERNAL = (
    :__has_mass_matrix, :_unwrap_val, :has_initializeprob, :parameterless_type,
)

# `@..` is owned by FastBroadcast and reaches StochasticDiffEqCore through
# DiffEqBase's re-export (the standard SciML access path); it is non-public in
# DiffEqBase, so it is ignored for the are-public check as well as via-owners.
const DIFFEQBASE_INTERNAL = (Symbol("@.."),)

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
