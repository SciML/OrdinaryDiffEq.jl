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

# Size dZ to carry MronRoe Fourier coefficients.
# MronRoe needs: 2*m*n + (m²+m)/2 random values.
# The noise process generates dZ as N(0,1) white noise of this length,
# which perform_step! unpacks into LevyAreaCoefficients.
function StochasticDiffEqCore._z_prototype(alg::RKMilGeneral, rand_prototype, iip::Bool)
    rand_prototype isa Number && return rand_prototype
    m = length(rand_prototype)
    # Determine truncation level
    p = alg.p
    if p === nothing
        # Default: use automatic truncation for a reasonable dt
        # This will be recomputed at solve time, but we need a size now.
        # Use a conservative n that covers typical step sizes.
        p = max(10, m)
    end
    n = p
    n_coeffs = StochasticDiffEqLevyArea.norv(m, n, MronRoe())
    rp2 = similar(rand_prototype, n_coeffs)
    rp2 .= false
    return rp2
end
