alg_order(alg::EM) = 1 // 2
alg_order(alg::LambaEM) = 1 // 2
alg_order(alg::SplitEM) = 1 // 2
alg_order(alg::PCEuler) = 1 // 2
alg_order(alg::EulerHeun) = 1 // 2
alg_order(alg::LambaEulerHeun) = 1 // 2
alg_order(alg::SimplifiedEM) = 1 // 2
alg_order(alg::RKMil) = 1 // 1
alg_order(alg::RKMilCommute) = 1 // 1

function SciMLBase.alg_interpretation(alg::EulerHeun)
    return SciMLBase.AlgorithmInterpretation.Stratonovich
end
function SciMLBase.alg_interpretation(alg::LambaEulerHeun)
    return SciMLBase.AlgorithmInterpretation.Stratonovich
end
function SciMLBase.alg_interpretation(alg::RKMil{interpretation}) where {interpretation}
    return interpretation
end
SciMLBase.alg_interpretation(alg::RKMilCommute) = alg.interpretation

alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::EM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::LambaEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::EulerHeun) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::LambaEulerHeun) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SplitEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::PCEuler) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::SimplifiedEM) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RKMil) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RKMilCommute) = true

is_split_step(::EM{split}) where {split} = split
is_split_step(::LambaEM{split}) where {split} = split

issplit(::SplitEM) = true
