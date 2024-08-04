## SciMLBase Trait Definitions
function SciMLBase.isautodifferentiable(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    true
end
function SciMLBase.allows_arbitrary_number_types(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    true
end
function SciMLBase.allowscomplex(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    true
end
function SciMLBase.forwarddiffs_model(alg::Union{OrdinaryDiffEqAdaptiveImplicitAlgorithm,
        DAEAlgorithm,
        OrdinaryDiffEqImplicitAlgorithm})
    alg_autodiff(alg) isa AutoForwardDiff
end
SciMLBase.forwarddiffs_model_time(alg::RosenbrockAlgorithm) = true

# isadaptive is defined below.

## OrdinaryDiffEq Internal Traits

isfsal(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = true
isfsal(tab::DiffEqBase.ExplicitRKTableau) = tab.fsal

# isfsal(alg::CompositeAlgorithm) = isfsal(alg.algs[alg.current])
# Pseudo Non-FSAL
#isfsal(alg::RKM) = false

isfirk(alg) = false

get_current_isfsal(alg, cache) = isfsal(alg)

dt_required(alg) = true

isdiscretealg(alg) = false

function alg_stability_size end
has_stiff_interpolation(alg) = false

# evaluates f(t[i])
_eval_index(f::F, t::Tuple{A}, _) where {F, A} = f(t[1])
function _eval_index(f::F, t::Tuple{A, Vararg}, i) where {F, A}
    if i == 1
        f(t[1])
    else
        _eval_index(f, Base.tail(t), i - 1)
    end
end

function get_current_isfsal(alg::CompositeAlgorithm, cache)
    _eval_index(isfsal, alg.algs, cache.current)::Bool
end

all_fsal(alg, cache) = isfsal(alg)
all_fsal(alg::CompositeAlgorithm, cache) = _all_fsal(alg.algs)

@generated function _all_fsal(algs::T) where {T <: Tuple}
    ex = Expr(:tuple, map(1:length(T.types)) do i
        :(isfsal(algs[$i]))
    end...)
    :(all($ex))
end

issplit(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false

function _composite_beta1_default(algs::Tuple{T1, T2}, current, ::Val{QT},
        beta2) where {T1, T2, QT}
    if current == 1
        return QT(beta1_default(algs[1], beta2))
    else
        return QT(beta1_default(algs[2], beta2))
    end
end

@generated function _composite_beta1_default(algs::T, current, ::Val{QT},
        beta2) where {T <: Tuple, QT}
    expr = Expr(:block)
    for i in 1:length(T.types)
        push!(expr.args, quote
            if current == $i
                return QT(beta1_default(algs[$i], beta2))
            end
        end)
    end
    return expr
end

function _composite_beta2_default(algs::Tuple{T1, T2}, current,
        ::Val{QT}) where {T1, T2, QT}
    if current == 1
        return QT(beta2_default(algs[1]))
    else
        return QT(beta2_default(algs[2]))
    end
end

@generated function _composite_beta2_default(algs::T, current,
        ::Val{QT}) where {T <: Tuple, QT}
    expr = Expr(:block)
    for i in 1:length(T.types)
        push!(expr.args, quote
            if current == $i
                return QT(beta2_default(algs[$i]))
            end
        end)
    end
    return expr
end

function fsal_typeof(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, rate_prototype)
    typeof(rate_prototype)
end

isimplicit(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
isimplicit(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm) = true
isimplicit(alg::OrdinaryDiffEqImplicitAlgorithm) = true
isimplicit(alg::CompositeAlgorithm) = any(isimplicit.(alg.algs))

isdtchangeable(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = true
isdtchangeable(alg::CompositeAlgorithm) = all(isdtchangeable.(alg.algs))

ismultistep(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
ismultistep(alg::CompositeAlgorithm) = any(ismultistep.(alg.algs))

isadaptive(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
isadaptive(alg::OrdinaryDiffEqAdaptiveAlgorithm) = true
isadaptive(alg::OrdinaryDiffEqCompositeAlgorithm) = all(isadaptive.(alg.algs))

has_special_newton_error(alg) = false

anyadaptive(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = isadaptive(alg)
anyadaptive(alg::OrdinaryDiffEqCompositeAlgorithm) = any(isadaptive, alg.algs)

has_dtnew_modification(alg) = false
dtnew_modification(integrator, alg, dtnew) = dtnew

isautoswitch(alg) = false
isautoswitch(alg::CompositeAlgorithm) = alg.choice_function isa AutoSwitch

function qmin_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    isadaptive(alg) ? 1 // 5 : 0
end
qmin_default(alg::CompositeAlgorithm) = maximum(qmin_default.(alg.algs))

qmax_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = 10
qmax_default(alg::CompositeAlgorithm) = minimum(qmax_default.(alg.algs))

function has_chunksize(alg::OrdinaryDiffEqAlgorithm)
    return alg isa Union{OrdinaryDiffEqExponentialAlgorithm,
        OrdinaryDiffEqAdaptiveExponentialAlgorithm,
        OrdinaryDiffEqImplicitAlgorithm,
        OrdinaryDiffEqAdaptiveImplicitAlgorithm,
        DAEAlgorithm,
        CompositeAlgorithm}
end
function get_chunksize(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have a chunk size defined.")
end
function get_chunksize(alg::Union{OrdinaryDiffEqExponentialAlgorithm{CS},
        OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS},
        OrdinaryDiffEqImplicitAlgorithm{CS},
        OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS},
        DAEAlgorithm{CS},
        CompositeAlgorithm{CS}}) where {CS}
    Val(CS)
end

function get_chunksize_int(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have a chunk size defined.")
end
function get_chunksize_int(alg::Union{OrdinaryDiffEqExponentialAlgorithm{CS},
        OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS},
        OrdinaryDiffEqImplicitAlgorithm{CS},
        OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS},
        DAEAlgorithm{CS},
        CompositeAlgorithm{CS}}) where {CS}
    CS
end
# get_chunksize(alg::CompositeAlgorithm) = get_chunksize(alg.algs[alg.current_alg])

function DiffEqBase.prepare_alg(
        alg::Union{
            OrdinaryDiffEqAdaptiveImplicitAlgorithm{0, AD,
                FDT},
            OrdinaryDiffEqImplicitAlgorithm{0, AD, FDT},
            DAEAlgorithm{0, AD, FDT},
            OrdinaryDiffEqExponentialAlgorithm{0, AD, FDT}},
        u0::AbstractArray{T},
        p, prob) where {AD, FDT, T}

    # If not using autodiff or norecompile mode or very large bitsize (like a dual number u0 already)
    # don't use a large chunksize as it will either error or not be beneficial
    if !(alg_autodiff(alg) isa AutoForwardDiff) ||
       (isbitstype(T) && sizeof(T) > 24) ||
       (prob.f isa ODEFunction &&
        prob.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper)
        return remake(alg, chunk_size = Val{1}())
    end

    L = StaticArrayInterface.known_length(typeof(u0))
    if L === nothing # dynamic sized
        # If chunksize is zero, pick chunksize right at the start of solve and
        # then do function barrier to infer the full solve
        x = if prob.f.colorvec === nothing
            length(u0)
        else
            maximum(prob.f.colorvec)
        end

        cs = ForwardDiff.pickchunksize(x)
        return remake(alg, chunk_size = Val{cs}())
    else # statically sized
        cs = pick_static_chunksize(Val{L}())
        return remake(alg, chunk_size = cs)
    end
end

# Linear Exponential doesn't have any of the AD stuff
function DiffEqBase.prepare_alg(
        alg::OrdinaryDiffEqLinearExponentialAlgorithm,
        u0::AbstractArray,
        p, prob)
    alg
end

@generated function pick_static_chunksize(::Val{chunksize}) where {chunksize}
    x = ForwardDiff.pickchunksize(chunksize)
    :(Val{$x}())
end

function DiffEqBase.prepare_alg(alg::CompositeAlgorithm, u0, p, prob)
    algs = map(alg -> DiffEqBase.prepare_alg(alg, u0, p, prob), alg.algs)
    CompositeAlgorithm(algs, alg.choice_function)
end

has_autodiff(alg::OrdinaryDiffEqAlgorithm) = false
function has_autodiff(alg::Union{
        OrdinaryDiffEqAdaptiveImplicitAlgorithm, OrdinaryDiffEqImplicitAlgorithm,
        CompositeAlgorithm, OrdinaryDiffEqExponentialAlgorithm, DAEAlgorithm})
    true
end

# Extract AD type parameter from algorithm, returning as Val to ensure type stability for boolean options.
function _alg_autodiff(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have an autodifferentiation option defined.")
end
_alg_autodiff(::OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD}) where {CS, AD} = Val{AD}()
_alg_autodiff(::DAEAlgorithm{CS, AD}) where {CS, AD} = Val{AD}()
_alg_autodiff(::OrdinaryDiffEqImplicitAlgorithm{CS, AD}) where {CS, AD} = Val{AD}()
_alg_autodiff(alg::CompositeAlgorithm) = _alg_autodiff(alg.algs[end])
function _alg_autodiff(::Union{OrdinaryDiffEqExponentialAlgorithm{CS, AD},
        OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD}
}) where {
        CS, AD
}
    Val{AD}()
end

function alg_autodiff(alg)
    autodiff = _alg_autodiff(alg)
    if autodiff == Val(false)
        return AutoFiniteDiff()
    elseif autodiff == Val(true)
        return AutoForwardDiff()
    else
        return _unwrap_val(autodiff)
    end
end

# end

# alg_autodiff(alg::CompositeAlgorithm) = alg_autodiff(alg.algs[alg.current_alg])
get_current_alg_autodiff(alg, cache) = alg_autodiff(alg)
function get_current_alg_autodiff(alg::CompositeAlgorithm, cache)
    _eval_index(alg_autodiff, alg.algs, cache.current)::Bool
end

function alg_difftype(alg::Union{
        OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD, FDT, ST, CJ
        },
        OrdinaryDiffEqImplicitAlgorithm{CS, AD, FDT, ST, CJ},
        OrdinaryDiffEqExponentialAlgorithm{CS, AD, FDT, ST, CJ},
        OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD, FDT, ST,
            CJ},
        DAEAlgorithm{CS, AD, FDT, ST, CJ}}) where {CS, AD, FDT, ST,
        CJ}
    FDT
end

function standardtag(alg::Union{
        OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD, FDT, ST, CJ
        },
        OrdinaryDiffEqImplicitAlgorithm{CS, AD, FDT, ST, CJ},
        OrdinaryDiffEqExponentialAlgorithm{CS, AD, FDT, ST, CJ},
        OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD, FDT, ST,
            CJ},
        DAEAlgorithm{CS, AD, FDT, ST, CJ}}) where {CS, AD, FDT, ST,
        CJ}
    ST
end

function concrete_jac(alg::Union{
        OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD, FDT, ST, CJ
        },
        OrdinaryDiffEqImplicitAlgorithm{CS, AD, FDT, ST, CJ},
        OrdinaryDiffEqExponentialAlgorithm{CS, AD, FDT, ST, CJ},
        OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD, FDT, ST,
            CJ},
        DAEAlgorithm{CS, AD, FDT, ST, CJ}}) where {CS, AD, FDT, ST,
        CJ}
    CJ
end

alg_extrapolates(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
alg_extrapolates(alg::CompositeAlgorithm) = any(alg_extrapolates.(alg.algs))
function alg_order(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    error("Order is not defined for this algorithm")
end
alg_order(alg::CompositeAlgorithm) = maximum(alg_order, alg.algs)

function get_current_alg_order(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, cache)
    alg_order(alg)
end
function get_current_alg_order(alg::CompositeAlgorithm, cache)
    _eval_index(alg_order, alg.algs, cache.current)::Int
end

get_current_alg_order(alg::OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm, cache) = cache.order
function get_current_adaptive_order(alg::OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm, cache)
    cache.order
end

#alg_adaptive_order(alg::OrdinaryDiffEqAdaptiveAlgorithm) = error("Algorithm is adaptive with no order")
function get_current_adaptive_order(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm},
        cache)
    alg_adaptive_order(alg)
end
function get_current_adaptive_order(alg::CompositeAlgorithm, cache)
    _eval_index(alg_adaptive_order, alg.algs, cache.current)::Int
end

alg_maximum_order(alg) = alg_order(alg)
alg_maximum_order(alg::CompositeAlgorithm) = maximum(alg_order(x) for x in alg.algs)

alg_adaptive_order(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = alg_order(alg) - 1

# this is actually incorrect and is purposefully decreased as this tends
# to track the real error much better
# this is actually incorrect and is purposefully decreased as this tends
# to track the real error much better

function default_controller(alg, cache, qoldinit, _beta1 = nothing, _beta2 = nothing)
    if ispredictive(alg)
        return PredictiveController()
    elseif isstandard(alg)
        return IController()
    else # Default is PI-controller
        QT = typeof(qoldinit)
        beta1, beta2 = _digest_beta1_beta2(alg, cache, Val(QT), _beta1, _beta2)
        return PIController(beta1, beta2)
    end
end

function _digest_beta1_beta2(alg, cache, ::Val{QT}, _beta1, _beta2) where {QT}
    if alg isa OrdinaryDiffEqCompositeAlgorithm
        beta2 = _beta2 === nothing ?
                _composite_beta2_default(alg.algs, cache.current, Val(QT)) : _beta2
        beta1 = _beta1 === nothing ?
                _composite_beta1_default(alg.algs, cache.current, Val(QT), beta2) : _beta1
    else
        beta2 = _beta2 === nothing ? beta2_default(alg) : _beta2
        beta1 = _beta1 === nothing ? beta1_default(alg, beta2) : _beta1
    end
    return convert(QT, beta1)::QT, convert(QT, beta2)::QT
end

# other special cases in controllers.jl
function beta2_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    isadaptive(alg) ? 2 // (5alg_order(alg)) : 0
end

function beta1_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, beta2)
    isadaptive(alg) ? 7 // (10alg_order(alg)) : 0
end

function gamma_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    isadaptive(alg) ? 9 // 10 : 0
end
gamma_default(alg::CompositeAlgorithm) = maximum(gamma_default, alg.algs)

fac_default_gamma(alg) = false

qsteady_min_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = 1
qsteady_max_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = 1
qsteady_max_default(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm) = 6 // 5
# But don't re-use Jacobian if not adaptive: too risky and cannot pull back
qsteady_max_default(alg::OrdinaryDiffEqImplicitAlgorithm) = isadaptive(alg) ? 1 // 1 : 0
#TODO
#DiffEqBase.nlsolve_default(::QNDF, ::Val{κ}) = 1//2

# SSP coefficients
ssp_coefficient(alg) = error("$alg is not a strong stability preserving method.")

# We shouldn't do this probably.
#ssp_coefficient(alg::ImplicitEuler) = Inf

alg_can_repeat_jac(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
alg_can_repeat_jac(alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm) = true

function unwrap_alg(alg::SciMLBase.DEAlgorithm, is_stiff)
    if !(alg isa CompositeAlgorithm)
        return alg
    elseif alg.choice_function isa AutoSwitchCache
        if length(alg.algs) > 2
            return alg.algs[alg.choice_function.current]
        end
        if is_stiff === nothing
            throwautoswitch(alg)
        end
        num = is_stiff ? 2 : 1
        if num == 1
            return alg.algs[1]
        else
            return alg.algs[2]
        end
    else
        error("this dispatch does not support this algorithm right now")
    end
end

function unwrap_alg(integrator, is_stiff)
    alg = integrator.alg
    if !(alg isa CompositeAlgorithm)
        return alg
    elseif alg.choice_function isa AutoSwitchCache
        if length(alg.algs) > 2
            alg.algs[alg.choice_function.current]
        else
            if is_stiff === nothing
                throwautoswitch(alg)
            end
            num = is_stiff ? 2 : 1
            if num == 1
                return alg.algs[1]
            else
                return alg.algs[2]
            end
        end
    else
        return _eval_index(identity, alg.algs, integrator.cache.current)
    end
end

function throwautoswitch(alg)
    throw(ArgumentError("one of $(alg.algs) is not compatible with stiffness-based autoswitching"))
end

# Whether `uprev` is used in the algorithm directly.
uses_uprev(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, adaptive::Bool) = true
uses_uprev(alg::OrdinaryDiffEqAdaptiveAlgorithm, adaptive::Bool) = true

ispredictive(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
ispredictive(alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm) = alg.controller === :Predictive
isstandard(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
isstandard(alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm) = alg.controller === :Standard

isWmethod(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false

isesdirk(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false

is_mass_matrix_alg(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
is_mass_matrix_alg(alg::CompositeAlgorithm) = all(is_mass_matrix_alg, alg.algs)
is_mass_matrix_alg(alg::RosenbrockAlgorithm) = true
is_mass_matrix_alg(alg::NewtonAlgorithm) = !isesdirk(alg)
