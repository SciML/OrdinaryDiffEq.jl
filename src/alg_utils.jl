## SciMLBase Trait Definitions

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
isfsal(alg::Rodas5) = false
isfsal(alg::Rodas5P) = false
isfsal(alg::Rodas4) = false
isfsal(alg::Rodas42) = false
isfsal(alg::Rodas4P) = false
isfsal(alg::Rodas4P2) = false
isfsal(alg::Vern7) = false
isfsal(alg::Vern8) = false
isfsal(alg::Vern9) = false
# Pseudo Non-FSAL
isfsal(alg::ORK256) = false
isfsal(alg::CarpenterKennedy2N54) = false
isfsal(alg::SHLDDRK64) = false
isfsal(alg::DGLDDRK73_C) = false
isfsal(alg::DGLDDRK84_C) = false
isfsal(alg::DGLDDRK84_F) = false
isfsal(alg::NDBLSRK124) = false
isfsal(alg::NDBLSRK134) = false
isfsal(alg::NDBLSRK144) = false
isfsal(alg::PDIRK44) = false
isfsal(alg::DImplicitEuler) = false
isfsal(alg::RKO65) = false
isfsal(alg::FRK65) = true
#isfsal(alg::RKM) = false

isfsal(alg::RDPK3Sp35) = false
isfsal(alg::RDPK3Sp49) = false
isfsal(alg::RDPK3Sp510) = false

isfsal(alg::SSPRK22) = false
isfsal(alg::SSPRK33) = false
isfsal(alg::SSPRK53) = false
isfsal(alg::SSPRK53_2N1) = false
isfsal(alg::SSPRK53_2N2) = false
isfsal(alg::SSPRK53_H) = false
isfsal(alg::SSPRK63) = false
isfsal(alg::SSPRK73) = false
isfsal(alg::SSPRK83) = false
isfsal(alg::SSPRK43) = false
isfsal(alg::SSPRK432) = false
isfsal(alg::SSPRK932) = false
isfsal(alg::SSPRK54) = false
isfsal(alg::SSPRK104) = false

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

function get_chunksize(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have a chunk size defined.")
end
get_chunksize(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD}) where {CS, AD} = Val(CS)
get_chunksize(alg::OrdinaryDiffEqImplicitAlgorithm{CS, AD}) where {CS, AD} = Val(CS)
get_chunksize(alg::DAEAlgorithm{CS, AD}) where {CS, AD} = Val(CS)
function get_chunksize(alg::Union{OrdinaryDiffEqExponentialAlgorithm{CS, AD},
    OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD}}) where {
    CS,
    AD,
}
    Val(CS)
end

function get_chunksize_int(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have a chunk size defined.")
end
get_chunksize_int(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD}) where {CS, AD} = CS
get_chunksize_int(alg::OrdinaryDiffEqImplicitAlgorithm{CS, AD}) where {CS, AD} = CS
get_chunksize_int(alg::DAEAlgorithm{CS, AD}) where {CS, AD} = CS
function get_chunksize_int(alg::Union{OrdinaryDiffEqExponentialAlgorithm{CS, AD},
    OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD}}) where {
    CS,
    AD,
}
    CS
end
# get_chunksize(alg::CompositeAlgorithm) = get_chunksize(alg.algs[alg.current_alg])

function DiffEqBase.prepare_alg(alg::Union{
        OrdinaryDiffEqAdaptiveImplicitAlgorithm{0, AD,
            FDT},
        OrdinaryDiffEqImplicitAlgorithm{0, AD, FDT},
        DAEAlgorithm{0, AD, FDT},
        OrdinaryDiffEqExponentialAlgorithm{0, AD, FDT}},
    u0::AbstractArray{T},
    p, prob) where {AD, FDT, T}
    if alg isa OrdinaryDiffEqExponentialAlgorithm
        linsolve = nothing
    elseif alg.linsolve === nothing
        if (prob.f isa ODEFunction && prob.f.f isa AbstractSciMLOperator)
            linsolve = LinearSolve.defaultalg(prob.f.f, u0)
        elseif (prob.f isa SplitFunction &&
                prob.f.f1.f isa AbstractSciMLOperator)
            linsolve = LinearSolve.defaultalg(prob.f.f1.f, u0)
            if (linsolve === nothing) || (linsolve isa LinearSolve.DefaultLinearSolver &&
                linsolve.alg !== LinearSolve.DefaultAlgorithmChoice.KrylovJL_GMRES)
                msg = "Split ODE problem do not work with factorization linear solvers. Bug detailed in https://github.com/SciML/OrdinaryDiffEq.jl/pull/1643. Defaulting to linsolve=KrylovJL()"
                @warn msg
                linsolve = KrylovJL_GMRES()
            end
        elseif (prob isa ODEProblem || prob isa DDEProblem) &&
               (prob.f.mass_matrix === nothing ||
                (prob.f.mass_matrix !== nothing &&
                 !(typeof(prob.f.jac_prototype) <: AbstractSciMLOperator)))
            linsolve = LinearSolve.defaultalg(prob.f.jac_prototype, u0)
        else
            # If mm is a sparse matrix and A is a MatrixOperator, then let linear
            # solver choose things later
            linsolve = nothing
        end
    else
        linsolve = alg.linsolve
    end

    # If norecompile mode or very large bitsize, like a dual number u0 already, then
    # don't use a large chunksize as it will either error or not be beneficial
    if (isbitstype(T) && sizeof(T) > 24) || (prob.f isa ODEFunction &&
        prob.f.f isa
        FunctionWrappersWrappers.FunctionWrappersWrapper)
        if alg isa OrdinaryDiffEqExponentialAlgorithm
            return remake(alg, chunk_size = Val{1}())
        else
            return remake(alg, chunk_size = Val{1}(), linsolve = linsolve)
        end
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

        if alg isa OrdinaryDiffEqExponentialAlgorithm
            return remake(alg, chunk_size = Val{cs}())
        else
            return remake(alg, chunk_size = Val{cs}(), linsolve = linsolve)
        end
    else # statically sized
        cs = pick_static_chunksize(Val{L}())
        if alg isa OrdinaryDiffEqExponentialAlgorithm
            return remake(alg, chunk_size = cs)
        else
            return remake(alg, chunk_size = cs, linsolve = linsolve)
        end
    end
end

# Linear Exponential doesn't have any of the AD stuff
function DiffEqBase.prepare_alg(alg::Union{ETD2, SplitEuler, LinearExponential,
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

# Extract AD type parameter from algorithm, returning as Val to ensure type stability for boolean options.
function _alg_autodiff(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have an autodifferentiation option defined.")
end
_alg_autodiff(::OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD}) where {CS, AD} = Val{AD}()
_alg_autodiff(::DAEAlgorithm{CS, AD}) where {CS, AD} = Val{AD}()
_alg_autodiff(::OrdinaryDiffEqImplicitAlgorithm{CS, AD}) where {CS, AD} = Val{AD}()
function _alg_autodiff(::Union{OrdinaryDiffEqExponentialAlgorithm{CS, AD},
    OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD},
}) where {
    CS, AD,
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
alg_extrapolates(alg::ImplicitEuler) = true
alg_extrapolates(alg::DImplicitEuler) = true
alg_extrapolates(alg::DABDF2) = true
alg_extrapolates(alg::Trapezoid) = true
alg_extrapolates(alg::SDIRK22) = true
alg_extrapolates(alg::IRKN4) = true
alg_extrapolates(alg::IRKN3) = true
alg_extrapolates(alg::ABDF2) = true
alg_extrapolates(alg::SBDF) = true
alg_extrapolates(alg::MEBDF2) = true
alg_extrapolates(alg::MagnusLeapfrog) = true

function alg_order(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    error("Order is not defined for this algorithm")
end
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
get_current_alg_order(alg::QNDF, cache) = cache.order
get_current_alg_order(alg::FBDF, cache) = cache.order
get_current_adaptive_order(alg::QNDF, cache) = cache.order
get_current_adaptive_order(alg::FBDF, cache) = cache.order
function get_current_adaptive_order(alg::OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm,
    cache)
    cache.cur_order
end
function get_current_alg_order(alg::OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm,
    cache)
    cache.cur_order
end
get_current_alg_order(alg::ExtrapolationMidpointDeuflhard, cache) = 2(cache.n_curr + 1)
get_current_alg_order(alg::ImplicitDeuflhardExtrapolation, cache) = 2(cache.n_curr + 1)
get_current_adaptive_order(alg::ExtrapolationMidpointDeuflhard, cache) = 2cache.n_curr
get_current_adaptive_order(alg::ImplicitDeuflhardExtrapolation, cache) = 2cache.n_curr
get_current_alg_order(alg::ExtrapolationMidpointHairerWanner, cache) = 2(cache.n_curr + 1)
get_current_alg_order(alg::ImplicitHairerWannerExtrapolation, cache) = 2(cache.n_curr + 1)
get_current_alg_order(alg::ImplicitEulerBarycentricExtrapolation, cache) = cache.n_curr
get_current_alg_order(alg::ImplicitEulerExtrapolation, cache) = cache.n_curr + 1
get_current_adaptive_order(alg::ExtrapolationMidpointHairerWanner, cache) = 2cache.n_curr
get_current_adaptive_order(alg::ImplicitHairerWannerExtrapolation, cache) = 2cache.n_curr
get_current_adaptive_order(alg::ImplicitEulerExtrapolation, cache) = cache.n_curr - 1
function get_current_adaptive_order(alg::ImplicitEulerBarycentricExtrapolation, cache)
    cache.n_curr - 2
end

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
alg_order(alg::AitkenNeville) = alg.init_order
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
alg_order(alg::RK46NL) = 4
alg_order(alg::KuttaPRK2p5) = 5
alg_order(alg::RKO65) = 5
alg_order(alg::FRK65) = 6

alg_order(alg::SymplecticEuler) = 1
alg_order(alg::VelocityVerlet) = 2
alg_order(alg::VerletLeapfrog) = 2
alg_order(alg::PseudoVerletLeapfrog) = 2
alg_order(alg::McAte2) = 2
alg_order(alg::Ruth3) = 3
alg_order(alg::McAte3) = 3
alg_order(alg::McAte4) = 4
alg_order(alg::CandyRoz4) = 4
alg_order(alg::CalvoSanz4) = 4
alg_order(alg::McAte42) = 4
alg_order(alg::McAte5) = 5
alg_order(alg::Yoshida6) = 6
alg_order(alg::KahanLi6) = 6
alg_order(alg::McAte8) = 8
alg_order(alg::KahanLi8) = 8
alg_order(alg::SofSpa10) = 10

alg_order(alg::IRKN3) = 3
alg_order(alg::Nystrom4) = 4
alg_order(alg::FineRKN5) = 5
alg_order(alg::Nystrom4VelocityIndependent) = 4
alg_order(alg::IRKN4) = 4
alg_order(alg::Nystrom5VelocityIndependent) = 5
alg_order(alg::DPRKN4) = 4
alg_order(alg::DPRKN5) = 5
alg_order(alg::DPRKN6) = 6
alg_order(alg::DPRKN6FM) = 6
alg_order(alg::DPRKN8) = 8
alg_order(alg::DPRKN12) = 12
alg_order(alg::ERKN4) = 4
alg_order(alg::ERKN5) = 5
alg_order(alg::ERKN7) = 7

alg_order(alg::Midpoint) = 2

alg_order(alg::ORK256) = 2
alg_order(alg::CarpenterKennedy2N54) = 4
alg_order(alg::SHLDDRK52) = 2
alg_order(alg::SHLDDRK_2N) = 4
alg_order(alg::SHLDDRK64) = 4
alg_order(alg::DGLDDRK73_C) = 3
alg_order(alg::DGLDDRK84_C) = 4
alg_order(alg::DGLDDRK84_F) = 4
alg_order(alg::NDBLSRK124) = 4
alg_order(alg::NDBLSRK134) = 4
alg_order(alg::NDBLSRK144) = 4
alg_order(alg::CFRLDDRK64) = 4
alg_order(alg::TSLDDRK74) = 4
alg_order(alg::CKLLSRK43_2) = 3
alg_order(alg::CKLLSRK54_3C) = 4
alg_order(alg::CKLLSRK95_4S) = 5
alg_order(alg::CKLLSRK95_4C) = 5
alg_order(alg::CKLLSRK95_4M) = 5
alg_order(alg::CKLLSRK54_3C_3R) = 4
alg_order(alg::CKLLSRK54_3M_3R) = 4
alg_order(alg::CKLLSRK54_3N_3R) = 4
alg_order(alg::CKLLSRK85_4C_3R) = 5
alg_order(alg::CKLLSRK85_4M_3R) = 5
alg_order(alg::CKLLSRK85_4P_3R) = 5
alg_order(alg::CKLLSRK54_3N_4R) = 4
alg_order(alg::CKLLSRK54_3M_4R) = 4
alg_order(alg::CKLLSRK65_4M_4R) = 5
alg_order(alg::CKLLSRK85_4FM_4R) = 5
alg_order(alg::CKLLSRK75_4M_5R) = 5
alg_order(alg::ParsaniKetchesonDeconinck3S32) = 2
alg_order(alg::ParsaniKetchesonDeconinck3S82) = 2
alg_order(alg::ParsaniKetchesonDeconinck3S53) = 3
alg_order(alg::ParsaniKetchesonDeconinck3S173) = 3
alg_order(alg::ParsaniKetchesonDeconinck3S94) = 4
alg_order(alg::ParsaniKetchesonDeconinck3S184) = 4
alg_order(alg::ParsaniKetchesonDeconinck3S105) = 5
alg_order(alg::ParsaniKetchesonDeconinck3S205) = 5
alg_order(alg::RDPK3Sp35) = 3
alg_order(alg::RDPK3SpFSAL35) = 3
alg_order(alg::RDPK3Sp49) = 4
alg_order(alg::RDPK3SpFSAL49) = 4
alg_order(alg::RDPK3Sp510) = 5
alg_order(alg::RDPK3SpFSAL510) = 5
alg_order(alg::KYK2014DGSSPRK_3S2) = 2

alg_order(alg::SSPRK22) = 2
alg_order(alg::SSPRKMSVS32) = 2
alg_order(alg::SSPRK33) = 3
alg_order(alg::KYKSSPRK42) = 2
alg_order(alg::SSPRK53) = 3
alg_order(alg::SSPRK53_2N1) = 3
alg_order(alg::SSPRK53_2N2) = 3
alg_order(alg::SSPRK53_H) = 3
alg_order(alg::SSPRK63) = 3
alg_order(alg::SSPRK73) = 3
alg_order(alg::SSPRK83) = 3
alg_order(alg::SSPRK43) = 3
alg_order(alg::SSPRK432) = 3
alg_order(alg::SSPRKMSVS43) = 3
alg_order(alg::SSPRK932) = 3
alg_order(alg::SSPRK54) = 4
alg_order(alg::SSPRK104) = 4

alg_order(alg::RK4) = 4
alg_order(alg::RKM) = 4
alg_order(alg::ExplicitRK) = alg.tableau.order
alg_order(alg::MSRK5) = 5
alg_order(alg::MSRK6) = 6
alg_order(alg::Stepanov5) = 5
alg_order(alg::SIR54) = 5

alg_order(alg::BS3) = 3
alg_order(alg::BS5) = 5
alg_order(alg::OwrenZen3) = 3
alg_order(alg::OwrenZen4) = 4
alg_order(alg::OwrenZen5) = 5
alg_order(alg::DP5) = 5
alg_order(alg::Tsit5) = 5
alg_order(alg::DP8) = 8
alg_order(alg::Vern6) = 6
alg_order(alg::Vern7) = 7
alg_order(alg::Vern8) = 8
alg_order(alg::Vern9) = 9
alg_order(alg::TanYam7) = 7
alg_order(alg::TsitPap8) = 8
alg_order(alg::RadauIIA3) = 3
alg_order(alg::RadauIIA5) = 5
alg_order(alg::ImplicitEuler) = 1
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
alg_order(alg::Trapezoid) = 2
alg_order(alg::ImplicitMidpoint) = 2
alg_order(alg::TRBDF2) = 2
alg_order(alg::SSPSDIRK2) = 2
alg_order(alg::SDIRK2) = 2
alg_order(alg::SDIRK22) = 2
alg_order(alg::Kvaerno3) = 3
alg_order(alg::Kvaerno4) = 4
alg_order(alg::Kvaerno5) = 5
alg_order(alg::ESDIRK54I8L2SA) = 5
alg_order(alg::ESDIRK436L2SA2) = 4
alg_order(alg::ESDIRK437L2SA) = 4
alg_order(alg::ESDIRK547L2SA2) = 5
alg_order(alg::ESDIRK659L2SA) = 6
alg_order(alg::KenCarp3) = 3
alg_order(alg::CFNLIRK3) = 3
alg_order(alg::KenCarp4) = 4
alg_order(alg::KenCarp47) = 4
alg_order(alg::KenCarp5) = 5
alg_order(alg::KenCarp58) = 5
alg_order(alg::Cash4) = 4
alg_order(alg::SFSDIRK4) = 4
alg_order(alg::SFSDIRK5) = 4
alg_order(alg::SFSDIRK6) = 4
alg_order(alg::SFSDIRK7) = 4
alg_order(alg::SFSDIRK8) = 4
alg_order(alg::Hairer4) = 4
alg_order(alg::Hairer42) = 4
alg_order(alg::Feagin10) = 10
alg_order(alg::Feagin12) = 12
alg_order(alg::Feagin14) = 14
alg_order(alg::PFRK87) = 8

alg_order(alg::Rosenbrock23) = 2
alg_order(alg::Rosenbrock32) = 3
alg_order(alg::ROS3P) = 3
alg_order(alg::Rodas3) = 3
alg_order(alg::ROS34PW1a) = 3
alg_order(alg::ROS34PW1b) = 3
alg_order(alg::ROS34PW2) = 3
alg_order(alg::ROS34PW3) = 4
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

alg_order(alg::ABDF2) = 2
alg_order(alg::QNDF1) = 1
alg_order(alg::QNDF2) = 2

alg_order(alg::QNDF) = 1 #dummy value
alg_order(alg::FBDF) = 1 #dummy value

alg_order(alg::SBDF) = alg.order

alg_order(alg::ROCK2) = 2
alg_order(alg::ROCK4) = 4

alg_order(alg::ESERK4) = 4
alg_order(alg::ESERK5) = 5
alg_order(alg::SERK2) = 2

alg_order(alg::RKC) = 2
alg_order(alg::IRKC) = 2

alg_order(alg::MEBDF2) = 2
alg_order(alg::PDIRK44) = 4

alg_order(alg::DImplicitEuler) = 1
alg_order(alg::DABDF2) = 2
alg_order(alg::DFBDF) = 1#dummy value

alg_order(alg::Alshina2) = 2
alg_order(alg::Alshina3) = 3
alg_order(alg::Alshina6) = 6

alg_maximum_order(alg) = alg_order(alg)
alg_maximum_order(alg::CompositeAlgorithm) = maximum(alg_order(x) for x in alg.algs)
alg_maximum_order(alg::ExtrapolationMidpointDeuflhard) = 2(alg.max_order + 1)
alg_maximum_order(alg::ImplicitDeuflhardExtrapolation) = 2(alg.max_order + 1)
alg_maximum_order(alg::ExtrapolationMidpointHairerWanner) = 2(alg.max_order + 1)
alg_maximum_order(alg::ImplicitHairerWannerExtrapolation) = 2(alg.max_order + 1)
alg_maximum_order(alg::ImplicitEulerExtrapolation) = 2(alg.max_order + 1)

alg_maximum_order(alg::ImplicitEulerBarycentricExtrapolation) = alg.max_order
alg_adaptive_order(alg::ExplicitRK) = alg.tableau.adaptiveorder
alg_adaptive_order(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = alg_order(alg) - 1
alg_adaptive_order(alg::Feagin10) = 8
alg_adaptive_order(alg::Feagin14) = 12

alg_adaptive_order(alg::Rosenbrock23) = 3
alg_adaptive_order(alg::Rosenbrock32) = 2

alg_adaptive_order(alg::RadauIIA3) = 1
alg_adaptive_order(alg::RadauIIA5) = 3
alg_adaptive_order(alg::RKC) = 2
alg_adaptive_order(alg::IRKC) = 1

alg_adaptive_order(alg::ImplicitEuler) = 0
alg_adaptive_order(alg::Trapezoid) = 1
# this is actually incorrect and is purposefully decreased as this tends
# to track the real error much better
alg_adaptive_order(alg::ImplicitMidpoint) = 1
# this is actually incorrect and is purposefully decreased as this tends
# to track the real error much better

alg_adaptive_order(alg::Exprb32) = 2
alg_adaptive_order(alg::Exprb43) = 4
alg_adaptive_order(alg::AN5) = 5

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

function default_controller(alg::Union{ExtrapolationMidpointDeuflhard,
        ImplicitDeuflhardExtrapolation,
        ExtrapolationMidpointHairerWanner,
        ImplicitHairerWannerExtrapolation,
        ImplicitEulerExtrapolation,
        ImplicitEulerBarycentricExtrapolation}, cache,
    qoldinit, _beta1 = nothing, _beta2 = nothing)
    QT = typeof(qoldinit)
    beta1, beta2 = _digest_beta1_beta2(alg, cache, Val(QT), _beta1, _beta2)
    return ExtrapolationController(beta1)
end

function _digest_beta1_beta2(alg, cache, ::Val{QT}, _beta1, _beta2) where {QT}
    if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
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

function default_controller(alg::RDPK3Sp35, cache, qoldinit, args...)
    QT = typeof(qoldinit)
    return PIDController(map(Base.Fix1(convert, QT), (0.64, -0.31, 0.04))...)
end

function default_controller(alg::RDPK3SpFSAL35, cache, qoldinit, args...)
    QT = typeof(qoldinit)
    return PIDController(map(Base.Fix1(convert, QT), (0.70, -0.23, 0.00))...)
end

function default_controller(alg::RDPK3Sp49, cache, qoldinit, args...)
    QT = typeof(qoldinit)
    return PIDController(map(Base.Fix1(convert, QT), (0.25, -0.12, 0.00))...)
end

function default_controller(alg::RDPK3SpFSAL49, cache, qoldinit, args...)
    QT = typeof(qoldinit)
    return PIDController(map(Base.Fix1(convert, QT), (0.38, -0.18, 0.01))...)
end

function default_controller(alg::RDPK3Sp510, cache, qoldinit, args...)
    QT = typeof(qoldinit)
    return PIDController(map(Base.Fix1(convert, QT), (0.47, -0.20, 0.06))...)
end

function default_controller(alg::RDPK3SpFSAL510, cache, qoldinit, args...)
    QT = typeof(qoldinit)
    return PIDController(map(Base.Fix1(convert, QT), (0.45, -0.13, 0.00))...)
end

# other special cases in controllers.jl
function default_controller(alg::Union{JVODE, QNDF, FBDF}, args...)
    DummyController()
end

function beta2_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    isadaptive(alg) ? 2 // (5alg_order(alg)) : 0
end
beta2_default(alg::FunctionMap) = 0
beta2_default(alg::DP8) = 0 // 1
beta2_default(alg::DP5) = 4 // 100
beta2_default(alg::ExtrapolationMidpointDeuflhard) = 0 // 1
beta2_default(alg::ImplicitDeuflhardExtrapolation) = 0 // 1
beta2_default(alg::ExtrapolationMidpointHairerWanner) = 0 // 1
beta2_default(alg::ImplicitHairerWannerExtrapolation) = 0 // 1
beta2_default(alg::ImplicitEulerExtrapolation) = 0 // 1
beta2_default(alg::ImplicitEulerBarycentricExtrapolation) = 0 // 1

function beta1_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, beta2)
    isadaptive(alg) ? 7 // (10alg_order(alg)) : 0
end
beta1_default(alg::FunctionMap, beta2) = 0
beta1_default(alg::DP8, beta2) = typeof(beta2)(1 // alg_order(alg)) - beta2 / 5
beta1_default(alg::DP5, beta2) = typeof(beta2)(1 // alg_order(alg)) - 3beta2 / 4
beta1_default(alg::ExtrapolationMidpointDeuflhard, beta2) = 1 // (2alg.init_order + 1)
beta1_default(alg::ImplicitDeuflhardExtrapolation, beta2) = 1 // (2alg.init_order + 1)
beta1_default(alg::ExtrapolationMidpointHairerWanner, beta2) = 1 // (2alg.init_order + 1)
beta1_default(alg::ImplicitHairerWannerExtrapolation, beta2) = 1 // (2alg.init_order + 1)
beta1_default(alg::ImplicitEulerExtrapolation, beta2) = 1 // (alg.init_order + 1)
beta1_default(alg::ImplicitEulerBarycentricExtrapolation, beta2) = 1 // (alg.init_order - 1)

function gamma_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    isadaptive(alg) ? 9 // 10 : 0
end
gamma_default(alg::RKC) = 8 // 10
gamma_default(alg::IRKC) = 8 // 10
function gamma_default(alg::ExtrapolationMidpointDeuflhard)
    (1 // 4)^beta1_default(alg, beta2_default(alg))
end
function gamma_default(alg::ImplicitDeuflhardExtrapolation)
    (1 // 4)^beta1_default(alg, beta2_default(alg))
end
function gamma_default(alg::ExtrapolationMidpointHairerWanner)
    (65 // 100)^beta1_default(alg, beta2_default(alg))
end
function gamma_default(alg::ImplicitHairerWannerExtrapolation)
    (65 // 100)^beta1_default(alg, beta2_default(alg))
end
function gamma_default(alg::ImplicitEulerExtrapolation)
    (65 // 100)^beta1_default(alg, beta2_default(alg))
end

function gamma_default(alg::ImplicitEulerBarycentricExtrapolation)
    (80 // 100)^beta1_default(alg, beta2_default(alg))
end

qsteady_min_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = 1
qsteady_min_default(alg::FBDF) = 9 // 10
qsteady_max_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = 1
qsteady_max_default(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm) = 6 // 5
# But don't re-use Jacobian if not adaptive: too risky and cannot pull back
qsteady_max_default(alg::OrdinaryDiffEqImplicitAlgorithm) = isadaptive(alg) ? 1 // 1 : 0
qsteady_max_default(alg::AN5) = 3 // 2
qsteady_max_default(alg::JVODE) = 3 // 2
qsteady_max_default(alg::QNDF1) = 2 // 1
qsteady_max_default(alg::QNDF2) = 2 // 1
qsteady_max_default(alg::QNDF) = 2 // 1
qsteady_max_default(alg::FBDF) = 2 // 1

#TODO
#DiffEqBase.nlsolve_default(::QNDF, ::Val{κ}) = 1//2

function FunctionMap_scale_by_time(alg::FunctionMap{scale_by_time}) where {scale_by_time}
    scale_by_time
end

# SSP coefficients
"""
    ssp_coefficient(alg)

Return the SSP coefficient of the ODE algorithm `alg`. If one time step of size
`dt` with `alg` can be written as a convex combination of explicit Euler steps
with step sizes `cᵢ * dt`, the SSP coefficient is the minimal value of `1/cᵢ`.

# Examples
```julia-repl
julia> ssp_coefficient(SSPRK104())
6
```
"""
ssp_coefficient(alg) = error("$alg is not a strong stability preserving method.")
ssp_coefficient(alg::Euler) = 1
ssp_coefficient(alg::SSPRK22) = 1
ssp_coefficient(alg::SSPRK33) = 1
ssp_coefficient(alg::KYKSSPRK42) = 2.459
ssp_coefficient(alg::SSPRK53) = 2.65
ssp_coefficient(alg::SSPRK53_2N1) = 2.18
ssp_coefficient(alg::SSPRK53_2N2) = 2.148
ssp_coefficient(alg::SSPRK53_H) = 2.65
ssp_coefficient(alg::SSPRK63) = 3.518
ssp_coefficient(alg::SSPRK73) = 4.2879
ssp_coefficient(alg::SSPRK83) = 5.107
ssp_coefficient(alg::SSPRK43) = 2
ssp_coefficient(alg::SSPRK432) = 2
ssp_coefficient(alg::SSPRKMSVS32) = 0.5
ssp_coefficient(alg::SSPRKMSVS43) = 0.33
ssp_coefficient(alg::SSPRK932) = 6
ssp_coefficient(alg::SSPRK54) = 1.508
ssp_coefficient(alg::SSPRK104) = 6
ssp_coefficient(alg::KYK2014DGSSPRK_3S2) = 0.8417

# We shouldn't do this probably.
#ssp_coefficient(alg::ImplicitEuler) = Inf
ssp_coefficient(alg::SSPSDIRK2) = 4

# stability regions
alg_stability_size(alg::ExplicitRK) = alg.tableau.stability_size
alg_stability_size(alg::DP5) = 3.3066
alg_stability_size(alg::Tsit5) = 3.5068
alg_stability_size(alg::Vern6) = 4.8553
alg_stability_size(alg::Vern7) = 4.6400
alg_stability_size(alg::Vern8) = 5.8641
alg_stability_size(alg::Vern9) = 4.4762

alg_can_repeat_jac(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
alg_can_repeat_jac(alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm) = true
alg_can_repeat_jac(alg::IRKC) = false

function unwrap_alg(alg::SciMLBase.DEAlgorithm, is_stiff)
    iscomp = typeof(alg) <: CompositeAlgorithm
    if !iscomp
        return alg
    elseif typeof(alg.choice_function) <: AutoSwitchCache
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
    iscomp = typeof(alg) <: CompositeAlgorithm
    if !iscomp
        return alg
    elseif typeof(alg.choice_function) <: AutoSwitchCache
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
        return _eval_index(identity, alg.algs, integrator.cache.current)
    end
end

function throwautoswitch(alg)
    throw(ArgumentError("one of $(alg.algs) is not compatible with stiffness-based autoswitching"))
end

# Whether `uprev` is used in the algorithm directly.
uses_uprev(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, adaptive::Bool) = true
uses_uprev(alg::ORK256, adaptive::Bool) = false
uses_uprev(alg::CarpenterKennedy2N54, adaptive::Bool) = false
uses_uprev(alg::SHLDDRK64, adaptive::Bool) = false
uses_uprev(alg::DGLDDRK73_C, adaptive::Bool) = false
uses_uprev(alg::DGLDDRK84_C, adaptive::Bool) = false
uses_uprev(alg::DGLDDRK84_F, adaptive::Bool) = false
uses_uprev(alg::NDBLSRK124, adaptive::Bool) = false
uses_uprev(alg::NDBLSRK134, adaptive::Bool) = false
uses_uprev(alg::NDBLSRK144, adaptive::Bool) = false
uses_uprev(alg::CFRLDDRK64, adaptive::Bool) = false
uses_uprev(alg::TSLDDRK74, adaptive::Bool) = false
uses_uprev(alg::OrdinaryDiffEqAdaptiveAlgorithm, adaptive::Bool) = true
uses_uprev(alg::CKLLSRK43_2, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3C, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK95_4S, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK95_4C, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK95_4M, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3C_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3M_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3N_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK85_4C_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK85_4M_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK85_4P_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3N_4R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3M_4R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK65_4M_4R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK85_4FM_4R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK75_4M_5R, adaptive::Bool) = adaptive

ispredictive(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
ispredictive(alg::Union{RKC}) = true
ispredictive(alg::Union{SERK2}) = alg.controller === :Predictive
ispredictive(alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm) = alg.controller === :Predictive
isstandard(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
isstandard(alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm) = alg.controller === :Standard
isstandard(alg::VCABM) = true

isWmethod(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
isWmethod(alg::Rosenbrock23) = true
isWmethod(alg::Rosenbrock32) = true
isWmethod(alg::ROS34PW1a) = true
isWmethod(alg::ROS34PW1b) = true
isWmethod(alg::ROS34PW2) = true
isWmethod(alg::ROS34PW3) = true
isWmethod(alg::RosenbrockW6S4OS) = true

isesdirk(alg::TRBDF2) = true
function isesdirk(alg::Union{KenCarp3, KenCarp4, KenCarp5, KenCarp58,
    Kvaerno3, Kvaerno4, Kvaerno5, ESDIRK437L2SA,
    ESDIRK54I8L2SA, ESDIRK436L2SA2, ESDIRK547L2SA2,
    ESDIRK659L2SA, CFNLIRK3})
    true
end
isesdirk(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false

is_mass_matrix_alg(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
is_mass_matrix_alg(alg::CompositeAlgorithm) = all(is_mass_matrix_alg, alg.algs)
is_mass_matrix_alg(alg::RosenbrockAlgorithm) = true
is_mass_matrix_alg(alg::NewtonAlgorithm) = !isesdirk(alg)
