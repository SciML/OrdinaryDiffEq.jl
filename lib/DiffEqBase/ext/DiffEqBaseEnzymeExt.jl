module DiffEqBaseEnzymeExt

@static if isempty(VERSION.prerelease)
    using DiffEqBase
    import DiffEqBase: value
    using Enzyme
    import Enzyme: Const, MixedDuplicated
    using ChainRulesCore


    @inline function copy_or_reuse(config, val, idx)
        if Enzyme.EnzymeRules.overwritten(config)[idx] && ismutable(val)
            return deepcopy(val)
        else
            return val
        end
    end

    @inline function arg_copy(data, i)
        config, args = data
        copy_or_reuse(config, args[i].val, i + 5)
    end

    # Note these following functions are generally not considered user facing from within Enzyme.
    # They enable additional performance/usability here (e.g. inactive kwargs).
    # Contact wsmoses@ before modifying (and beware their semantics may change without semver).

    Enzyme.EnzymeRules.inactive_kwarg(::typeof(DiffEqBase.solve_up), prob, sensalg::Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm}, u0, p, args...; kwargs...) = nothing

    Enzyme.EnzymeRules.has_easy_rule(::typeof(DiffEqBase.solve_up), prob, sensalg::Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm}, u0, p, args...; kwargs...) = nothing

    function Enzyme.EnzymeRules.augmented_primal(
            config::Enzyme.EnzymeRules.RevConfigWidth{1},
            func::Const{typeof(DiffEqBase.solve_up)}, RTA::Union{Type{Duplicated{RT}}, Type{MixedDuplicated{RT}}}, prob,
            sensealg::Union{
                Const{Nothing}, Const{<:DiffEqBase.AbstractSensitivityAlgorithm},
            },
            u0, p, args...; kwargs...
        ) where {RT}

        res = DiffEqBase._solve_adjoint(
            copy_or_reuse(config, prob.val, 2), copy_or_reuse(config, sensealg.val, 3),
            copy_or_reuse(config, u0.val, 4), copy_or_reuse(config, p.val, 5),
            SciMLBase.EnzymeOriginator(), ntuple(Base.Fix1(arg_copy, (config, args)), Val(length(args)))...;
            kwargs...
        )

        primal = if Enzyme.EnzymeRules.needs_primal(config)
            res[1]
        else
            nothing
        end

        shadow = if Enzyme.EnzymeRules.needs_shadow(config)
            Enzyme.make_zero(res[1])::RT
        else
            nothing
        end
        tup = if Enzyme.EnzymeRules.needs_shadow(config)
            (shadow, res[2])
        else
            nothing
        end
        return Enzyme.EnzymeRules.augmented_rule_return_type(config, RTA)(primal, shadow, tup)
    end

    function Enzyme.EnzymeRules.reverse(
            config::Enzyme.EnzymeRules.RevConfigWidth{1},
            func::Const{typeof(DiffEqBase.solve_up)}, ::Union{Type{Duplicated{RT}}, Type{MixedDuplicated{RT}}}, tape, prob,
            sensealg::Union{
                Const{Nothing}, Const{<:DiffEqBase.AbstractSensitivityAlgorithm},
            },
            u0, p, args...; kwargs...
        ) where {RT}

        if Enzyme.EnzymeRules.needs_shadow(config)
            dres, clos = tape
            dres = dres::RT
            dargs = clos(dres)
            for (darg, ptr) in zip(dargs, (func, prob, sensealg, u0, p, args...))
                if ptr isa Enzyme.Const
                    continue
                end
                if darg == ChainRulesCore.NoTangent()
                    continue
                end
                if ptr isa MixedDuplicated
                    ptr.dval[] .+= darg
                else
                    ptr.dval .+= darg
                end
            end
            Enzyme.make_zero!(dres.u)
        end

        return ntuple(Returns(nothing), Val(length(args) + 4))
    end
end

end
