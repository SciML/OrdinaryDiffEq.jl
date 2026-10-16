using ExplicitImports
using OrdinaryDiffEqCore
using Reactant
using Test

const ReactantExt = Base.get_extension(OrdinaryDiffEqCore, :OrdinaryDiffEqCoreReactantExt)

@test ReactantExt isa Module

# Reactant has no public equivalent of these two names (see the extension); any other
# non-public Reactant name fails here.
@test check_all_explicit_imports_are_public(
    ReactantExt; ignore = (:TracedRNumber, :bitcast_convert)
) === nothing

@test check_all_qualified_accesses_are_public(
    ReactantExt; ignore = (:controller_fastpower, :dt_below_time_eps, :value_eps)
) === nothing
