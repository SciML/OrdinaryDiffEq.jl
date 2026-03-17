alg_order(alg::ImplicitEM) = 1 // 2
alg_order(alg::ImplicitEulerHeun) = 1 // 2
alg_order(alg::ImplicitRKMil) = 1 // 1
alg_order(alg::ISSEM) = 1 // 2
alg_order(alg::ISSEulerHeun) = 1 // 2
alg_order(alg::SKenCarp) = 2 // 1

alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::ImplicitEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::ImplicitEulerHeun) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::ISSEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::ISSEulerHeun) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SKenCarp) = true
function alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::ImplicitRKMil)
    return is_diagonal_noise(prob)
end

alg_needs_extra_process(alg::SKenCarp) = true

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

function SciMLBase.alg_interpretation(alg::ImplicitEulerHeun)
    return SciMLBase.AlgorithmInterpretation.Stratonovich
end
