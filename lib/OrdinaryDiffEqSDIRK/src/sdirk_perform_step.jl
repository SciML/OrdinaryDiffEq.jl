# Legacy SDIRK perform_step implementations
# 
# These have been replaced by the unified tableau-based generic implementations in:
# - generic_sdirk_perform_step.jl (contains the unified perform_step! functions)
# - unified_sdirk_tableaus.jl (contains the unified tableau structures)
#
# All SDIRK methods now use the generic implementation with their respective tableaus,
# eliminating code duplication while maintaining full backward compatibility and performance.