## alg_order

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

alg_order(alg::IRI1) = 1 // 1  # Strong order, weak order is 2

## alg_interpretation (Stratonovich methods)

SciMLBase.alg_interpretation(alg::RS1) = SciMLBase.AlgorithmInterpretation.Stratonovich
SciMLBase.alg_interpretation(alg::RS2) = SciMLBase.AlgorithmInterpretation.Stratonovich

SciMLBase.alg_interpretation(alg::NON) = SciMLBase.AlgorithmInterpretation.Stratonovich
SciMLBase.alg_interpretation(alg::COM) = SciMLBase.AlgorithmInterpretation.Stratonovich
SciMLBase.alg_interpretation(alg::NON2) = SciMLBase.AlgorithmInterpretation.Stratonovich

## alg_compatible

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
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::IRI1) = true

## alg_needs_extra_process

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

## _z_prototype overrides (Z process sizing for algorithms with extra noise)

# PL1WM needs m*(m-1)/2 auxiliary random variables for iterated integrals
function _z_prototype(alg::PL1WM, rand_prototype, iip::Bool)
    if !iip && rand_prototype isa Number
        return nothing
    end
    m = length(rand_prototype)
    rp2 = similar(rand_prototype, Int(m * (m - 1) / 2))
    rp2 .= false
    return rp2
end

# W2Ito1 always needs exactly 2 auxiliary random variables
function _z_prototype(alg::W2Ito1, rand_prototype, iip::Bool)
    return zeros(eltype(rand_prototype), 2)
end
