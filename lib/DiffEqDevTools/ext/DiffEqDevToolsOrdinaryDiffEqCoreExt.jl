module DiffEqDevToolsOrdinaryDiffEqCoreExt

import DiffEqDevTools
import OrdinaryDiffEqCore

const ODEC = OrdinaryDiffEqCore

function DiffEqDevTools._preset_algorithm_tags(
        alg::Union{
            ODEC.OrdinaryDiffEqAlgorithm,
            ODEC.DAEAlgorithm,
            ODEC.StochasticDiffEqAlgorithm,
            ODEC.StochasticDiffEqRODEAlgorithm,
        }
    )
    tags = Symbol[]

    if alg isa ODEC.DAEAlgorithm
        append!(tags, (:dae, :implicit))
    elseif alg isa ODEC.OrdinaryDiffEqAlgorithm
        push!(tags, ODEC.isimplicit(alg) ? :implicit : :explicit)
    elseif alg isa ODEC.StochasticDiffEqAlgorithm
        push!(tags, :sde)
    else
        push!(tags, :rode)
    end

    alg isa ODEC.RosenbrockAlgorithm && push!(tags, :rosenbrock)
    alg isa ODEC.ExponentialAlgorithm && push!(tags, :exponential)
    alg isa ODEC.PartitionedAlgorithm && push!(tags, :partitioned)
    alg isa ODEC.OrdinaryDiffEqCompositeAlgorithm && push!(tags, :composite)
    if alg isa ODEC.OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm
        append!(tags, (:adams, :variable_order, :multistep))
    end

    applicable(ODEC.isfirk, alg) && ODEC.isfirk(alg) && push!(tags, :firk)
    if applicable(ODEC.isesdirk, alg) && ODEC.isesdirk(alg)
        append!(tags, (:sdirk, :esdirk))
    end
    applicable(ODEC.issplit, alg) && ODEC.issplit(alg) && push!(tags, :split)
    applicable(ODEC.ismultistep, alg) && ODEC.ismultistep(alg) &&
        push!(tags, :multistep)

    return unique!(tags)
end

end
