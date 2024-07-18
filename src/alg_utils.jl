## SciMLBase Trait Definitions
using OrdinaryDiffEq: AbstractController
function SciMLBase.isautodifferentiable(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm,
        FunctionMap})
    true
end
function SciMLBase.allows_arbitrary_number_types(alg::Union{OrdinaryDiffEqAlgorithm,
        DAEAlgorithm, FunctionMap})
    true
end
function SciMLBase.allowscomplex(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm,
        FunctionMap})
    true
end
SciMLBase.isdiscrete(alg::FunctionMap) = true
function SciMLBase.forwarddiffs_model(alg::Union{OrdinaryDiffEqAdaptiveImplicitAlgorithm,
        DAEAlgorithm,
        OrdinaryDiffEqImplicitAlgorithm,
        ExponentialAlgorithm})
    alg_autodiff(alg) isa AutoForwardDiff
end
SciMLBase.forwarddiffs_model_time(alg::RosenbrockAlgorithm) = true

# isadaptive is defined below.

## OrdinaryDiffEq Internal Traits

isfsal(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = true
isfsal(tab::DiffEqBase.ExplicitRKTableau) = tab.fsal

# isfsal(alg::CompositeAlgorithm) = isfsal(alg.algs[alg.current])
isfsal(alg::FunctionMap) = false
isfsal(alg::Rodas3P) = false
isfsal(alg::Rodas23W) = false
isfsal(alg::Rodas5) = false
isfsal(alg::Rodas5P) = false
isfsal(alg::Rodas5Pr) = false
isfsal(alg::Rodas5Pe) = false
isfsal(alg::Rodas4) = false
isfsal(alg::Rodas42) = false
isfsal(alg::Rodas4P) = false
isfsal(alg::Rodas4P2) = false
# Pseudo Non-FSAL
isfsal(alg::PDIRK44) = false
isfsal(alg::DImplicitEuler) = false
isfsal(alg::RKO65) = false
isfsal(alg::FRK65) = true
#isfsal(alg::RKM) = false

isfsal(alg::PSRK3p5q4) = false
isfsal(alg::PSRK3p6q5) = false
isfsal(alg::PSRK4p7q6) = false

get_current_isfsal(alg, cache) = isfsal(alg)

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
issplit(alg::SplitAlgorithms) = true

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
fsal_typeof(alg::ETD2, rate_prototype) = ETD2Fsal{typeof(rate_prototype)}
function fsal_typeof(alg::CompositeAlgorithm, rate_prototype)
    fsal = map(x -> fsal_typeof(x, rate_prototype), alg.algs)
    @assert length(unique(fsal))==1 "`fsal_typeof` must be consistent"
    return fsal[1]
end

isimplicit(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
isimplicit(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm) = true
isimplicit(alg::OrdinaryDiffEqImplicitAlgorithm) = true
isimplicit(alg::CompositeAlgorithm) = any(isimplicit.(alg.algs))

isdtchangeable(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = true
isdtchangeable(alg::CompositeAlgorithm) = all(isdtchangeable.(alg.algs))

function isdtchangeable(alg::Union{LawsonEuler, NorsettEuler, LieEuler, MagnusGauss4,
        CayleyEuler, ETDRK2, ETDRK3, ETDRK4, HochOst4, ETD2})
    false
end # due to caching

ismultistep(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
ismultistep(alg::CompositeAlgorithm) = any(ismultistep.(alg.algs))
ismultistep(alg::ETD2) = true

isadaptive(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
isadaptive(alg::OrdinaryDiffEqAdaptiveAlgorithm) = true
isadaptive(alg::OrdinaryDiffEqCompositeAlgorithm) = all(isadaptive.(alg.algs))
isadaptive(alg::DImplicitEuler) = true
isadaptive(alg::DABDF2) = true
isadaptive(alg::DFBDF) = true

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
qmin_default(alg::DP8) = 1 // 3

qmax_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = 10
qmax_default(alg::CompositeAlgorithm) = minimum(qmax_default.(alg.algs))
qmax_default(alg::DP8) = 6
qmax_default(alg::Union{RadauIIA3, RadauIIA5}) = 8

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
        alg::Union{ETD2, SplitEuler, LinearExponential,
            OrdinaryDiffEqLinearExponentialAlgorithm},
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
alg_extrapolates(alg::DImplicitEuler) = true
alg_extrapolates(alg::DABDF2) = true
alg_extrapolates(alg::MagnusLeapfrog) = true

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
get_current_alg_order(alg::JVODE, cache) = get_current_adaptive_order(alg, cache)

#alg_adaptive_order(alg::OrdinaryDiffEqAdaptiveAlgorithm) = error("Algorithm is adaptive with no order")
function get_current_adaptive_order(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm},
        cache)
    alg_adaptive_order(alg)
end
function get_current_adaptive_order(alg::CompositeAlgorithm, cache)
    _eval_index(alg_adaptive_order, alg.algs, cache.current)::Int
end

alg_order(alg::FunctionMap) = 0
alg_order(alg::Euler) = 1
alg_order(alg::Heun) = 2
alg_order(alg::Ralston) = 2
alg_order(alg::LawsonEuler) = 1
alg_order(alg::NorsettEuler) = 1
alg_order(alg::LieEuler) = 1
alg_order(alg::CayleyEuler) = 2
alg_order(alg::ETDRK2) = 2
alg_order(alg::ETDRK3) = 3
alg_order(alg::ETDRK4) = 4
alg_order(alg::HochOst4) = 4
alg_order(alg::Exp4) = 4
alg_order(alg::EPIRK4s3A) = 4
alg_order(alg::EPIRK4s3B) = 4
alg_order(alg::EPIRK5s3) = 5
alg_order(alg::EPIRK5P1) = 5
alg_order(alg::EPIRK5P2) = 5
alg_order(alg::EXPRB53s3) = 5
alg_order(alg::SplitEuler) = 1
alg_order(alg::ETD2) = 2
alg_order(alg::Exprb32) = 3
alg_order(alg::Exprb43) = 4
alg_order(alg::Anas5) = 5
alg_order(alg::KuttaPRK2p5) = 5
alg_order(alg::RKO65) = 5
alg_order(alg::FRK65) = 6

alg_order(alg::Midpoint) = 2

alg_order(alg::RK4) = 4
alg_order(alg::RKM) = 4
alg_order(alg::ExplicitRK) = alg.tableau.order
alg_order(alg::MSRK5) = 5
alg_order(alg::MSRK6) = 6
alg_order(alg::Stepanov5) = 5
alg_order(alg::SIR54) = 5
alg_order(alg::PSRK4p7q6) = 4
alg_order(alg::PSRK3p6q5) = 3
alg_order(alg::PSRK3p5q4) = 3

alg_order(alg::BS3) = 3
alg_order(alg::BS5) = 5
alg_order(alg::OwrenZen3) = 3
alg_order(alg::OwrenZen4) = 4
alg_order(alg::OwrenZen5) = 5
alg_order(alg::DP5) = 5
alg_order(alg::Tsit5) = 5
alg_order(alg::DP8) = 8
alg_order(alg::TanYam7) = 7
alg_order(alg::TsitPap8) = 8
alg_order(alg::RadauIIA3) = 3
alg_order(alg::RadauIIA5) = 5
alg_order(alg::RKMK2) = 2
alg_order(alg::RKMK4) = 4
alg_order(alg::LieRK4) = 4
alg_order(alg::CG2) = 2
alg_order(alg::CG3) = 3
alg_order(alg::CG4a) = 4
alg_order(alg::MagnusMidpoint) = 2
alg_order(alg::MagnusGauss4) = 4
alg_order(alg::MagnusNC6) = 6
alg_order(alg::MagnusGL6) = 6
alg_order(alg::MagnusGL8) = 8
alg_order(alg::MagnusNC8) = 8
alg_order(alg::MagnusGL4) = 4
alg_order(alg::MagnusAdapt4) = 4
alg_order(alg::LinearExponential) = 1
alg_order(alg::MagnusLeapfrog) = 2
alg_order(alg::PFRK87) = 8

alg_order(alg::ROS2) = 2
alg_order(alg::ROS2PR) = 2
alg_order(alg::ROS2S) = 2
alg_order(alg::ROS3) = 3
alg_order(alg::ROS3PR) = 3
alg_order(alg::Scholz4_7) = 3
alg_order(alg::Rosenbrock23) = 2
alg_order(alg::Rodas23W) = 3
alg_order(alg::Rosenbrock32) = 3
alg_order(alg::ROS3P) = 3
alg_order(alg::Rodas3) = 3
alg_order(alg::Rodas3P) = 3
alg_order(alg::ROS34PW1a) = 3
alg_order(alg::ROS34PW1b) = 3
alg_order(alg::ROS34PW2) = 3
alg_order(alg::ROS34PW3) = 4
alg_order(alg::ROS34PRw) = 3
alg_order(alg::ROS3PRL) = 3
alg_order(alg::ROS3PRL2) = 3
alg_order(alg::ROK4a) = 4
alg_order(alg::RosShamp4) = 4
alg_order(alg::Veldd4) = 4
alg_order(alg::Velds4) = 4
alg_order(alg::GRK4T) = 4
alg_order(alg::GRK4A) = 4
alg_order(alg::Ros4LStab) = 4
alg_order(alg::RosenbrockW6S4OS) = 4
alg_order(alg::Rodas4) = 4
alg_order(alg::Rodas42) = 4
alg_order(alg::Rodas4P) = 4
alg_order(alg::Rodas4P2) = 4
alg_order(alg::Rodas5) = 5
alg_order(alg::Rodas5P) = 5
alg_order(alg::Rodas5Pr) = 5
alg_order(alg::Rodas5Pe) = 5

alg_order(alg::AB3) = 3
alg_order(alg::AB4) = 4
alg_order(alg::AB5) = 5
alg_order(alg::ABM32) = 3
alg_order(alg::ABM43) = 4
alg_order(alg::ABM54) = 5

alg_order(alg::VCAB3) = 3
alg_order(alg::VCAB4) = 4
alg_order(alg::VCAB5) = 5
alg_order(alg::VCABM3) = 3
alg_order(alg::VCABM4) = 4
alg_order(alg::VCABM5) = 5

alg_order(alg::VCABM) = 1  #dummy value

alg_order(alg::CNAB2) = 2
alg_order(alg::CNLF2) = 2

alg_order(alg::AN5) = 5
alg_order(alg::JVODE) = 1  #dummy value


alg_order(alg::PDIRK44) = 4

alg_order(alg::DImplicitEuler) = 1
alg_order(alg::DABDF2) = 2
alg_order(alg::DFBDF) = 1#dummy value

alg_order(alg::Alshina2) = 2
alg_order(alg::Alshina3) = 3
alg_order(alg::Alshina6) = 6

alg_order(alg::QPRK98) = 9

alg_maximum_order(alg) = alg_order(alg)
alg_maximum_order(alg::CompositeAlgorithm) = maximum(alg_order(x) for x in alg.algs)

alg_adaptive_order(alg::ExplicitRK) = alg.tableau.adaptiveorder
alg_adaptive_order(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = alg_order(alg) - 1

alg_adaptive_order(alg::Rosenbrock23) = 3
alg_adaptive_order(alg::Rosenbrock32) = 2

alg_adaptive_order(alg::RadauIIA3) = 1
alg_adaptive_order(alg::RadauIIA5) = 3

# this is actually incorrect and is purposefully decreased as this tends
# to track the real error much better
# this is actually incorrect and is purposefully decreased as this tends
# to track the real error much better

alg_adaptive_order(alg::Exprb32) = 2
alg_adaptive_order(alg::Exprb43) = 4
alg_adaptive_order(alg::AN5) = 5

struct DummyController <: AbstractController
end

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
function default_controller(alg::Union{JVODE}, args...)
    DummyController()
end

function beta2_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    isadaptive(alg) ? 2 // (5alg_order(alg)) : 0
end
beta2_default(alg::FunctionMap) = 0
beta2_default(alg::DP8) = 0 // 1
beta2_default(alg::DP5) = 4 // 100

function beta1_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, beta2)
    isadaptive(alg) ? 7 // (10alg_order(alg)) : 0
end
beta1_default(alg::FunctionMap, beta2) = 0
beta1_default(alg::DP8, beta2) = typeof(beta2)(1 // alg_order(alg)) - beta2 / 5
beta1_default(alg::DP5, beta2) = typeof(beta2)(1 // alg_order(alg)) - 3beta2 / 4

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
qsteady_max_default(alg::AN5) = 3 // 2
qsteady_max_default(alg::JVODE) = 3 // 2

#TODO
#DiffEqBase.nlsolve_default(::QNDF, ::Val{κ}) = 1//2

function FunctionMap_scale_by_time(alg::FunctionMap{scale_by_time}) where {scale_by_time}
    scale_by_time
end

# SSP coefficients
ssp_coefficient(alg) = error("$alg is not a strong stability preserving method.")
ssp_coefficient(alg::Euler) = 1

# We shouldn't do this probably.
#ssp_coefficient(alg::ImplicitEuler) = Inf

# stability regions
alg_stability_size(alg::ExplicitRK) = alg.tableau.stability_size
alg_stability_size(alg::DP5) = 3.3066
alg_stability_size(alg::Tsit5) = 3.5068

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
isstandard(alg::VCABM) = true

isWmethod(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
isWmethod(alg::Rosenbrock23) = true
isWmethod(alg::Rosenbrock32) = true
isWmethod(alg::Rodas23W) = true
isWmethod(alg::ROS2S) = true
isWmethod(alg::ROS34PW1a) = true
isWmethod(alg::ROS34PW1b) = true
isWmethod(alg::ROS34PW2) = true
isWmethod(alg::ROS34PW3) = true
isWmethod(alg::ROS34PRw) = true
isWmethod(alg::ROK4a) = true
isWmethod(alg::RosenbrockW6S4OS) = true

isesdirk(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false

is_mass_matrix_alg(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
is_mass_matrix_alg(alg::CompositeAlgorithm) = all(is_mass_matrix_alg, alg.algs)
is_mass_matrix_alg(alg::RosenbrockAlgorithm) = true
is_mass_matrix_alg(alg::NewtonAlgorithm) = !isesdirk(alg)