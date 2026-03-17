alg_order(alg::RKMilGeneral) = 1 // 1
alg_order(alg::WangLi3SMil_A) = 1 // 1
alg_order(alg::WangLi3SMil_B) = 1 // 1
alg_order(alg::WangLi3SMil_C) = 1 // 1
alg_order(alg::WangLi3SMil_D) = 1 // 1
alg_order(alg::WangLi3SMil_E) = 1 // 1
alg_order(alg::WangLi3SMil_F) = 1 // 1

SciMLBase.alg_interpretation(alg::RKMilGeneral) = alg.interpretation

alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::RKMilGeneral) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_A) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_B) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_C) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_D) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_E) = true
alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::WangLi3SMil_F) = true

alg_needs_extra_process(alg::RKMilGeneral) = true

# Override Z prototype for RKMilGeneral — needs extra Z when p is specified
function StochasticDiffEqCore._z_prototype(alg::RKMilGeneral, rand_prototype, iip::Bool)
    if alg.p === nothing || rand_prototype isa Number
        return nothing
    end
    m = length(rand_prototype)
    rp2 = similar(rand_prototype, Int(m + alg.p * m * 2))
    rp2 .= false
    return rp2
end
