alg_order(alg::IIF1M) = 1 // 2
alg_order(alg::IIF2M) = 1 // 2
alg_order(alg::IIF1Mil) = 1 // 1

alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::IIF1M) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::IIF2M) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::IIF1Mil) = true
