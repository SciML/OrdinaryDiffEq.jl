## SciMLBase Trait Definitions

SciMLBase.allows_late_binding_tstops(::StochasticDiffEqAlgorithm) = true
SciMLBase.allows_late_binding_tstops(::StochasticDiffEqRODEAlgorithm) = true

SciMLBase.supports_solve_rng(
    ::SciMLBase.AbstractSDEProblem,
    ::StochasticDiffEqAlgorithm,
) = true

SciMLBase.supports_solve_rng(
    ::SciMLBase.AbstractRODEProblem,
    ::StochasticDiffEqRODEAlgorithm,
) = true

function SciMLBase.isautodifferentiable(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
            StochasticDiffEqJumpAlgorithm,
        }
    )
    return true
end
function SciMLBase.allows_arbitrary_number_types(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
            StochasticDiffEqJumpAlgorithm,
        }
    )
    return true
end
function SciMLBase.allowscomplex(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
            StochasticDiffEqJumpAlgorithm,
        }
    )
    return true
end
SciMLBase.isdiscrete(alg::StochasticDiffEqJumpAlgorithm) = true

function SciMLBase.forwarddiffs_model(
        alg::Union{
            StochasticDiffEqNewtonAlgorithm,
            StochasticDiffEqNewtonAdaptiveAlgorithm, StochasticDiffEqJumpNewtonAdaptiveAlgorithm,
            StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm,
        }
    )
    return OrdinaryDiffEqCore.alg_autodiff(alg) isa ADTypes.AutoForwardDiff
end

# Required for initialization, because ODECore._initialize_dae! calls it during
# OverrideInit
function OrdinaryDiffEqCore.has_autodiff(
        ::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
            StochasticDiffEqJumpAlgorithm,
        }
    )
    return false
end
for T in [
        StochasticDiffEqNewtonAlgorithm, StochasticDiffEqNewtonAdaptiveAlgorithm,
        StochasticDiffEqJumpNewtonAdaptiveAlgorithm,
        StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm,
    ]
    @eval OrdinaryDiffEqCore.has_autodiff(::$T) = true
end

_alg_autodiff(::StochasticDiffEqNewtonAlgorithm{T, AD}) where {T, AD} = Val{AD}()
_alg_autodiff(::StochasticDiffEqNewtonAdaptiveAlgorithm{T, AD}) where {T, AD} = Val{AD}()
function _alg_autodiff(::StochasticDiffEqJumpNewtonAdaptiveAlgorithm{T, AD}) where {T, AD}
    return Val{AD}()
end
function _alg_autodiff(
        ::StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm{
            T, AD,
        }
    ) where {T, AD}
    return Val{AD}()
end
_alg_autodiff(alg::StochasticCompositeAlgorithm) = _alg_autodiff(alg.algs[end])

function OrdinaryDiffEqCore.alg_autodiff(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    ad = _alg_autodiff(alg)
    if ad == Val(false)
        return ADTypes.AutoFiniteDiff()
    elseif ad == Val(true)
        return ADTypes.AutoForwardDiff()
    else
        return SciMLBase._unwrap_val(ad)
    end
end

isadaptive(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = false
function isadaptive(
        alg::Union{
            StochasticDiffEqAdaptiveAlgorithm, StochasticDiffEqRODEAdaptiveAlgorithm,
            StochasticDiffEqJumpAdaptiveAlgorithm,
            StochasticDiffEqJumpDiffusionAdaptiveAlgorithm,
        }
    )
    return true
end
function isadaptive(
        alg::Union{
            StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
        }
    )
    return all(isadaptive.(alg.algs))
end
function isadaptive(
        prob, alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    return isadaptive(alg)
end
isadaptive(prob::JumpProblem, alg::ImplicitEM) = false

## StochasticDiffEq Internal Traits

function qmax_default(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm})
    return isadaptive(alg) ? 9 // 8 : 0
end
function qmin_default(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm})
    return isadaptive(alg) ? 1 // 5 : 0
end

delta_default(alg) = 1 // 1
delta_default(alg::SRIW1) = 1 // 6

ispredictive(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = false
isstandard(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = false
function qsteady_min_default(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    return 1
end
function qsteady_max_default(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    return 1
end

# special cases in stepsize_controllers.jl
function default_controller(
        alg::Union{TauLeaping, CaoTauLeaping, ThetaTrapezoidalTauLeaping}, args...
    )
    return DummyController()
end

# Extend ODE's isaposteriori for SDE algorithms.
# CaoTauLeaping always accepts (a posteriori dt), others use error estimates.
OrdinaryDiffEqCore.isaposteriori(alg::CaoTauLeaping) = true

# SDE algorithms never extrapolate (no uprev2 storage).
OrdinaryDiffEqCore.alg_extrapolates(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = false

# SDE algorithms are never FSAL (update_fsal! is a no-op for SDE).
OrdinaryDiffEqCore.isfsal(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = false

alg_order(alg::EM) = 1 // 2
alg_order(alg::LambaEM) = 1 // 2
alg_order(alg::ImplicitEM) = 1 // 2
alg_order(alg::ImplicitEulerHeun) = 1 // 2
alg_order(alg::ImplicitRKMil) = 1 // 1
alg_order(alg::WangLi3SMil_A) = 1 // 1
alg_order(alg::WangLi3SMil_B) = 1 // 1
alg_order(alg::WangLi3SMil_C) = 1 // 1
alg_order(alg::WangLi3SMil_D) = 1 // 1
alg_order(alg::WangLi3SMil_E) = 1 // 1
alg_order(alg::WangLi3SMil_F) = 1 // 1
alg_order(alg::ISSEM) = 1 // 2
alg_order(alg::ISSEulerHeun) = 1 // 2
alg_order(alg::SplitEM) = 1 // 2
alg_order(alg::PCEuler) = 1 // 2
alg_order(alg::IIF1M) = 1 // 2
alg_order(alg::IIF2M) = 1 // 2
alg_order(alg::IIF1Mil) = 1 // 1
alg_order(alg::EulerHeun) = 1 // 2
alg_order(alg::LambaEulerHeun) = 1 // 2
alg_order(alg::RandomEM) = 1 // 2
alg_order(alg::RandomHeun) = 1 // 2
alg_order(alg::RandomTamedEM) = 1 // 2
alg_order(alg::SimplifiedEM) = 1 // 2
alg_order(alg::RKMil) = 1 // 1
alg_order(alg::RKMilCommute) = 1 // 1
alg_order(alg::RKMilGeneral) = 1 // 1

# Generalised version of SROCK1, both Ito ans Stratonovich, will have strong order of 1//2
# and weak order of 1 for Multidimensional Weiner process
# Stratonovich version strong order 1 for 1 dimensional Weiner Process or if noise is commutative
# Ito version can have strong order version for 1 dimensional Weiner Process,
# diagonal noise or commutative noise
alg_order(alg::SROCK1) = 1 // 2
alg_order(alg::SROCK2) = 1 // 1
alg_order(alg::KomBurSROCK2) = 1 // 1
alg_order(alg::SROCKC2) = 1 // 1
alg_order(alg::SROCKEM) = alg.strong_order_1 ? 1 // 1 : 1 // 2
alg_order(alg::SKSROCK) = 1 // 2
alg_order(alg::TangXiaoSROCK2) = 1 // 1

alg_order(alg::SRI) = alg.tableau.order
alg_order(alg::SRIW1) = 3 // 2
alg_order(alg::SRIW2) = 3 // 2
alg_order(alg::SOSRI) = 3 // 2
alg_order(alg::SOSRI2) = 3 // 2
alg_order(alg::SRA) = alg.tableau.order
alg_order(alg::SRA1) = 2 // 1
alg_order(alg::SRA2) = 2 // 1
alg_order(alg::SRA3) = 2 // 1
alg_order(alg::SOSRA) = 2 // 1
alg_order(alg::SOSRA2) = 2 // 1

alg_order(alg::DRI1) = 1 // 1
alg_order(alg::DRI1NM) = 1 // 1
alg_order(alg::RI1) = 1 // 1
alg_order(alg::RI3) = 1 // 1
alg_order(alg::RI5) = 1 // 1
alg_order(alg::RI6) = 1 // 1
alg_order(alg::RDI1WM) = 1 // 1
alg_order(alg::RDI2WM) = 1 // 1
alg_order(alg::RDI3WM) = 1 // 1
alg_order(alg::RDI4WM) = 1 // 1
alg_order(alg::W2Ito1) = 1 // 1

alg_order(alg::RS1) = 1 // 1
alg_order(alg::RS2) = 1 // 1

alg_order(alg::PL1WM) = 1 // 1
alg_order(alg::PL1WMA) = 1 // 1

alg_order(alg::NON) = 1 // 1
alg_order(alg::COM) = 1 // 1
alg_order(alg::NON2) = 1 // 1

alg_order(alg::SIEA) = 1 // 1
alg_order(alg::SMEA) = 1 // 1
alg_order(alg::SIEB) = 1 // 1
alg_order(alg::SMEB) = 1 // 1

alg_order(alg::TauLeaping) = 1 // 1
alg_order(alg::CaoTauLeaping) = 1 // 1
alg_order(alg::ImplicitTauLeaping) = 1 // 1  # Weak first order (backward Euler)
alg_order(alg::ThetaTrapezoidalTauLeaping) = 2 // 1  # Weak second order

alg_order(alg::BAOAB) = 1 // 1

alg_order(alg::SKenCarp) = 2 // 1
alg_order(alg::IRI1) = 1 // 1  # Strong order, weak order is 2
function alg_order(
        alg::Union{
            StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
        }
    )
    return maximum(alg_order.(alg.algs))
end
get_current_alg_order(alg::StochasticDiffEqAlgorithm, cache) = alg_order(alg)
function get_current_alg_order(
        alg::Union{
            StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
        },
        cache
    )
    return alg_order(alg.algs[cache.current])
end

function beta2_default(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm})
    return isadaptive(alg) ? 2 // (5alg_order(alg)) : 0
end
function beta1_default(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}, beta2)
    return isadaptive(alg) ? 7 // (10alg_order(alg)) : 0
end

isdtchangeable(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = true

function SciMLBase.alg_interpretation(alg::StochasticDiffEqAlgorithm)
    return SciMLBase.AlgorithmInterpretation.Ito
end
function SciMLBase.alg_interpretation(alg::EulerHeun)
    return SciMLBase.AlgorithmInterpretation.Stratonovich
end
function SciMLBase.alg_interpretation(alg::LambaEulerHeun)
    return SciMLBase.AlgorithmInterpretation.Stratonovich
end
function SciMLBase.alg_interpretation(alg::KomBurSROCK2)
    return SciMLBase.AlgorithmInterpretation.Stratonovich
end
function SciMLBase.alg_interpretation(alg::RKMil{interpretation}) where {interpretation}
    return interpretation
end
function SciMLBase.alg_interpretation(
        alg::SROCK1{
            interpretation, E,
        }
    ) where {interpretation, E}
    return interpretation
end
SciMLBase.alg_interpretation(alg::RKMilCommute) = alg.interpretation
SciMLBase.alg_interpretation(alg::RKMilGeneral) = alg.interpretation
function SciMLBase.alg_interpretation(
        alg::ImplicitRKMil{
            CS,
            AD,
            F,
            P,
            FDT,
            ST,
            CJ,
            N,
            T2,
            Controller,
            interpretation,
        }
    ) where {
        CS, AD, F, P, FDT, ST, CJ, N, T2, Controller, interpretation,
    }
    return interpretation
end

SciMLBase.alg_interpretation(alg::RS1) = SciMLBase.AlgorithmInterpretation.Stratonovich
SciMLBase.alg_interpretation(alg::RS2) = SciMLBase.AlgorithmInterpretation.Stratonovich

SciMLBase.alg_interpretation(alg::NON) = SciMLBase.AlgorithmInterpretation.Stratonovich
SciMLBase.alg_interpretation(alg::COM) = SciMLBase.AlgorithmInterpretation.Stratonovich
SciMLBase.alg_interpretation(alg::NON2) = SciMLBase.AlgorithmInterpretation.Stratonovich

function alg_compatible(
        prob, alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    return true
end
alg_compatible(prob, alg::StochasticDiffEqAlgorithm) = false

function alg_compatible(prob::JumpProblem, alg::StochasticDiffEqAlgorithm)
    return alg_compatible(prob.prob, alg) && prob.regular_jump === nothing &&
        prob.prob isa DiffEqBase.AbstractSDEProblem
end
alg_compatible(prob::JumpProblem, alg::EM) = alg_compatible(prob.prob, alg)
alg_compatible(prob::JumpProblem, alg::ImplicitEM) = alg_compatible(prob.prob, alg)

alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SRI) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SRIW1) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SRIW2) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SOSRI) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SOSRI2) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SRA) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SRA1) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SRA2) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SRA3) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SOSRA) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SOSRA2) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::DRI1) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::DRI1NM) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RI1) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RI3) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RI5) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RI6) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RDI1WM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RDI2WM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RDI3WM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RDI4WM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::W2Ito1) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RS1) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RS2) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::PL1WM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::PL1WMA) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::NON) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::COM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::NON2) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SIEA) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SMEA) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SIEB) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SMEB) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SKenCarp) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::IRI1) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::EM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::LambaEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_A) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_B) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_C) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_D) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_E) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_F) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SROCK1) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SROCK2) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::KomBurSROCK2) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SROCKC2) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SROCKEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SKSROCK) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::TangXiaoSROCK2) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::EulerHeun) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::LambaEulerHeun) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SplitEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::PCEuler) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::ImplicitEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::ImplicitEulerHeun) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::ISSEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::ISSEulerHeun) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SimplifiedEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RKMil) = is_diagonal_noise(prob)
function alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::ImplicitRKMil)
    return is_diagonal_noise(prob)
end
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RKMilCommute) = true # No good check for commutative noise
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RKMilGeneral) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::IIF1M) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::IIF2M) = true
function alg_compatible(
        prob::DiffEqBase.AbstractSDEProblem,
        alg::Union{
            StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
        }
    )
    return max((alg_compatible(prob, a) for a in alg.algs)...)
end
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::BAOAB) = is_diagonal_noise(prob)

function alg_compatible(
        prob::JumpProblem,
        alg::Union{StochasticDiffEqJumpAdaptiveAlgorithm, StochasticDiffEqJumpAlgorithm}
    )
    return prob.prob isa DiscreteProblem
end

function alg_needs_extra_process(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    return false
end
function alg_needs_extra_process(
        alg::Union{
            StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
        }
    )
    return max((alg_needs_extra_process(a) for a in alg.algs)...)
end
alg_needs_extra_process(alg::RKMilGeneral) = true
alg_needs_extra_process(alg::SRI) = true
alg_needs_extra_process(alg::SRIW1) = true
alg_needs_extra_process(alg::SRIW2) = true
alg_needs_extra_process(alg::SOSRI) = true
alg_needs_extra_process(alg::SOSRI2) = true
alg_needs_extra_process(alg::SRA) = true
alg_needs_extra_process(alg::SRA1) = true
alg_needs_extra_process(alg::SRA2) = true
alg_needs_extra_process(alg::SRA3) = true
alg_needs_extra_process(alg::SOSRA) = true
alg_needs_extra_process(alg::SOSRA2) = true
alg_needs_extra_process(alg::SKenCarp) = true
alg_needs_extra_process(alg::IRI1) = true
alg_needs_extra_process(alg::DRI1) = true
alg_needs_extra_process(alg::RI1) = true
alg_needs_extra_process(alg::RI3) = true
alg_needs_extra_process(alg::RI5) = true
alg_needs_extra_process(alg::RI6) = true
alg_needs_extra_process(alg::RDI1WM) = true
alg_needs_extra_process(alg::RDI2WM) = true
alg_needs_extra_process(alg::RDI3WM) = true
alg_needs_extra_process(alg::RDI4WM) = true
alg_needs_extra_process(alg::W2Ito1) = true
alg_needs_extra_process(alg::RS1) = true
alg_needs_extra_process(alg::RS2) = true
alg_needs_extra_process(alg::PL1WM) = true
alg_needs_extra_process(alg::NON) = true
alg_needs_extra_process(alg::NON2) = true

function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqNewtonAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end
function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqNewtonAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end
function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqJumpNewtonAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end
function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end

function OrdinaryDiffEqCore.get_current_alg_autodiff(alg::StochasticDiffEqCompositeAlgorithm, cache)
    return OrdinaryDiffEqCore.alg_autodiff(alg.algs[cache.current])
end

function OrdinaryDiffEqCore.get_chunksize(
        alg::StochasticDiffEqNewtonAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val(CS)
end
function OrdinaryDiffEqCore.get_chunksize(
        alg::StochasticDiffEqNewtonAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val(CS)
end
function OrdinaryDiffEqCore.get_chunksize(
        alg::StochasticDiffEqJumpNewtonAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val(CS)
end
function OrdinaryDiffEqCore.get_chunksize(
        alg::StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val(CS)
end

@static if isdefined(OrdinaryDiffEqCore, :standardtag)
    OrdinaryDiffEqCore.standardtag(
        alg::Union{
            StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller},
            StochasticDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ, Controller},
        }
    ) where {CS, AD, FDT, ST, CJ, Controller} = ST
end

@static if isdefined(OrdinaryDiffEqCore, :alg_difftype)
    OrdinaryDiffEqCore.alg_difftype(
        alg::Union{
            StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller},
            StochasticDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ, Controller},
        }
    ) where {
        CS, AD, FDT, ST, CJ, Controller,
    } = FDT
end

@static if isdefined(OrdinaryDiffEqCore, :concrete_jac)
    OrdinaryDiffEqCore.concrete_jac(
        alg::Union{
            StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller},
            StochasticDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ, Controller},
        }
    ) where {
        CS, AD, FDT, ST, CJ, Controller,
    } = CJ
end

alg_mass_matrix_compatible(alg::StochasticDiffEqAlgorithm) = false
alg_mass_matrix_compatible(alg::StochasticDiffEqRODEAlgorithm) = false
alg_can_repeat_jac(alg::StochasticDiffEqAlgorithm) = true

function alg_mass_matrix_compatible(
        alg::Union{
            StochasticDiffEqNewtonAlgorithm, StochasticDiffEqNewtonAdaptiveAlgorithm,
        }
    )
    if alg.symplectic
        return true
    elseif alg.theta == 1
        return true
    else
        error("Algorithm must be set as symplectic or theta=1 for mass matrices")
    end
end

is_split_step(::StochasticDiffEqAlgorithm) = false
is_split_step(::EM{split}) where {split} = split
is_split_step(::LambaEM{split}) where {split} = split

alg_stability_size(alg::SOSRI2) = 10.6
alg_stability_size(alg::SOSRA2) = 5.3

is_composite(alg) = false
is_composite(alg::StochasticDiffEqCompositeAlgorithm) = true
is_composite(alg::StochasticDiffEqRODECompositeAlgorithm) = true
function unwrap_alg(integrator, is_nlsolve)
    alg = integrator.alg
    if !is_composite(alg)
        return alg
    elseif alg.choice_function isa AutoSwitchCache
        num = is_nlsolve ? 2 : 1
        if num == 1
            return alg.algs[1]
        elseif num == 2
            return alg.algs[2]
        else
            return alg.algs[num]
        end
    else
        if integrator.cache.current == 1
            return alg.algs[1]
        elseif integrator.cache.current == 2
            return alg.algs[2]
        else
            return alg.algs[integrator.cache.current]
        end
    end
end

issplit(::StochasticDiffEqAlgorithm) = false
issplit(::SplitSDEAlgorithms) = true

function OrdinaryDiffEqCore.unwrap_alg(integrator::SDEIntegrator, is_stiff)
    alg = integrator.alg
    if !is_composite(alg)
        return alg
    elseif alg.choice_function isa AutoSwitchCache
        num = is_stiff ? 2 : 1
        if num == 1
            return alg.algs[1]
        elseif num == 2
            return alg.algs[2]
        else
            return alg.algs[num]
        end
    else
        if integrator.cache.current == 1
            return alg.algs[1]
        elseif integrator.cache.current == 2
            return alg.algs[2]
        else
            return alg.algs[integrator.cache.current]
        end
    end
end

alg_control_rate(::StochasticDiffEqAlgorithm) = false
alg_control_rate(::StochasticDiffEqRODEAlgorithm) = false
alg_control_rate(::TauLeaping) = true
alg_control_rate(::ThetaTrapezoidalTauLeaping) = true  # We manage rates manually
