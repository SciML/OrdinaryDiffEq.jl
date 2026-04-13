## alg_order
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

## delta_default
delta_default(alg::SRIW1) = 1 // 6

## alg_compatible
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

## alg_needs_extra_process
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

## alg_stability_size
alg_stability_size(alg::SOSRI2) = 10.6
alg_stability_size(alg::SOSRA2) = 5.3
