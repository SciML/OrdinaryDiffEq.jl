# Legacy SDIRK perform_step implementations
# 
# This file previously contained ~3100 lines of duplicated SDIRK perform_step! implementations.
# These have been replaced by the unified tableau-based generic implementations in:
# - generic_sdirk_perform_step.jl (contains the unified perform_step! functions)
# - unified_sdirk_tableaus.jl (contains the unified tableau structures)
#
# All SDIRK methods now use the generic implementation with their respective tableaus,
# eliminating code duplication while maintaining full backward compatibility and performance.
#
# The old implementations that were removed included:
# - ImplicitEuler, ImplicitMidpoint, Trapezoidal (basic methods)
# - SDIRK2, SDIRK22, SSPSDIRK2, Cash4 (2nd order methods)  
# - SFSDIRK4/5/6/7/8 (higher order methods)
# - ESDIRK54I8L2SA, ESDIRK436L2SA2, ESDIRK437L2SA, ESDIRK547L2SA2, ESDIRK659L2SA (embedded methods)
# - Hairer4 (Hairer's method)
# - All KenCarp and Kvaerno methods with additive splitting
#
# Special cases like TRBDF2 maintain their unique implementations in generic_sdirk_perform_step.jl
# but use the unified tableau structure.

# Re-export the generic implementations
# (These are automatically available since generic_sdirk_perform_step.jl is included)