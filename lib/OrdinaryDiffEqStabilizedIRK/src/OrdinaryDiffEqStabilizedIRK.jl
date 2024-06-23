module OrdinaryDiffEqStabilizedRK

import OrdinaryDiffEq: alg_order, alg_maximum_order,
                       calculate_residuals!,                       
                       beta2_default, beta1_default, gamma_default, issplit,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals,
                       OrdinaryDiffEqAlgorithm, OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqAdaptiveImplicitAlgorithm,
                       alg_cache, _unwrap_val,
                        _reshape, _vec, NLNewton, update_W!,
                       build_nlsolver, markfirststage!, du_alias_or_new,
                       nlsolve!, isnewton

using DiffEqBase, FastBroadcast, MuladdMacro, RecursiveArrayTools
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA

macro cache(expr)
    name = expr.args[2].args[1].args[1]
    fields = [x for x in expr.args[3].args if typeof(x) != LineNumberNode]
    cache_vars = Expr[]
    jac_vars = Pair{Symbol, Expr}[]
    for x in fields
        if x.args[2] == :uType || x.args[2] == :rateType ||
           x.args[2] == :kType || x.args[2] == :uNoUnitsType
            push!(cache_vars, :(c.$(x.args[1])))
        elseif x.args[2] == :DiffCacheType
            push!(cache_vars, :(c.$(x.args[1]).du))
            push!(cache_vars, :(c.$(x.args[1]).dual_du))
        end
    end
    quote
        $(esc(expr))
        $(esc(:full_cache))(c::$name) = tuple($(cache_vars...))
    end
end

include("algorithms.jl")
include("alg_utils.jl")
include("irkc_utils.jl")
include("irkc_caches.jl")
include("irkc_perform_step.jl")

export IRKC

end
