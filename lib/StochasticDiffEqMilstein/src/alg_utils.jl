alg_order(alg::RKMilGeneral) = 1 // 1
alg_order(alg::WangLi3SMil_A) = 1 // 1
alg_order(alg::WangLi3SMil_B) = 1 // 1
alg_order(alg::WangLi3SMil_C) = 1 // 1
alg_order(alg::WangLi3SMil_D) = 1 // 1
alg_order(alg::WangLi3SMil_E) = 1 // 1
alg_order(alg::WangLi3SMil_F) = 1 // 1

SciMLBase.alg_interpretation(alg::RKMilGeneral) = alg.interpretation

alg_compatible(prob::SciMLBase.AbstractSDEProblem, alg::RKMilGeneral) = true
alg_compatible(prob::SciMLBase.AbstractSDEProblem, alg::WangLi3SMil_A) = true
alg_compatible(prob::SciMLBase.AbstractSDEProblem, alg::WangLi3SMil_B) = true
alg_compatible(prob::SciMLBase.AbstractSDEProblem, alg::WangLi3SMil_C) = true
alg_compatible(prob::SciMLBase.AbstractSDEProblem, alg::WangLi3SMil_D) = true
alg_compatible(prob::SciMLBase.AbstractSDEProblem, alg::WangLi3SMil_E) = true
alg_compatible(prob::SciMLBase.AbstractSDEProblem, alg::WangLi3SMil_F) = true

alg_needs_extra_process(alg::RKMilGeneral) = true

"""
    _levyarea_terms(alg::RKMilGeneral, m, dt) -> n

Number of Fourier terms to retain in the Lévy area series for `m`-dimensional noise at
step size `dt`, honouring an explicit `alg.p` when the user set one.

`perform_step!` always evaluates the Lévy area with `MronRoe()`, which is order 1 in
`n`, so the count must come from that algorithm's error bound and the same accuracy
target `ε = c dt^(γ + 1/2)` that `get_iterated_I!` uses. That takes
`O(dt^(-1/2))` terms, not the `O(dt^(-1))` of an order-1/2 Fourier expansion. Getting
this wrong is expensive in memory and not just in `randn` calls, because the
coefficients are drawn as the noise process' `dZ`, so every step of them is kept in the
saved noise history.

Falls back to a fixed count when the step size is not a real number known up front —
the Z prototype has to be built before an automatically chosen initial step size is
known, and unitful step sizes do not go through `terms_needed`.
"""
function _levyarea_terms(alg::RKMilGeneral, m, dt)
    alg.p isa Integer && return max(1, alg.p)
    (dt isa Real && isfinite(dt) && !iszero(dt)) || return max(10, m)
    h = abs(float(dt))
    ε = alg.c * h^(alg_order(alg) + 1 // 2)
    n = StochasticDiffEqLevyArea.terms_needed(
        m, h, ε, MronRoe(), StochasticDiffEqLevyArea.MaxL2()
    )
    return max(1, n)
end

# Size dZ to carry MronRoe Fourier coefficients.
# MronRoe needs: 2*m*n + (m²+m)/2 random values.
# The noise process generates dZ as N(0,1) white noise of this length,
# which perform_step! unpacks into LevyAreaCoefficients.
function StochasticDiffEqCore._z_prototype(alg::RKMilGeneral, rand_prototype, iip::Bool)
    return StochasticDiffEqCore._z_prototype(alg, rand_prototype, iip, nothing)
end

function StochasticDiffEqCore._z_prototype(alg::RKMilGeneral, rand_prototype, iip::Bool, dt)
    rand_prototype isa Number && return rand_prototype
    m = length(rand_prototype)
    n_coeffs = StochasticDiffEqLevyArea.norv(m, _levyarea_terms(alg, m, dt), MronRoe())
    rp2 = similar(rand_prototype, n_coeffs)
    rp2 .= false
    return rp2
end
