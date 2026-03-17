# Generalised version of SROCK1, both Ito and Stratonovich, will have strong order of 1//2
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

function SciMLBase.alg_interpretation(alg::KomBurSROCK2)
    return SciMLBase.AlgorithmInterpretation.Stratonovich
end
function SciMLBase.alg_interpretation(
        alg::SROCK1{
            interpretation, E,
        }
    ) where {interpretation, E}
    return interpretation
end

alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SROCK1) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SROCK2) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::KomBurSROCK2) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SROCKC2) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SROCKEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SKSROCK) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::TangXiaoSROCK2) = true
