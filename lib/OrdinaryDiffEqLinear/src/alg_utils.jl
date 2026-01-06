alg_order(alg::MagnusMidpoint) = 2
alg_order(alg::RKMK2) = 2
alg_order(alg::RKMK4) = 4
alg_order(alg::LieRK4) = 4
alg_order(alg::CG3) = 3
alg_order(alg::CG2) = 2
alg_order(alg::CG4a) = 4
alg_order(alg::MagnusAdapt4) = 4
alg_order(alg::MagnusGauss4) = 4
alg_order(alg::MagnusNC6) = 6
alg_order(alg::MagnusGL6) = 6
alg_order(alg::MagnusGL8) = 8
alg_order(alg::MagnusNC8) = 8
alg_order(alg::MagnusGL4) = 4
alg_order(alg::LinearExponential) = 1
alg_order(alg::MagnusLeapfrog) = 2
alg_order(alg::LieEuler) = 1
alg_order(alg::CayleyEuler) = 2

alg_extrapolates(alg::MagnusLeapfrog) = true

dt_required(alg::LinearExponential) = false

function DiffEqBase.prepare_alg(
        alg::LinearExponential,
        u0::AbstractArray,
        p, prob
    )
    return alg
end

function isdtchangeable(alg::Union{LieEuler, MagnusGauss4, CayleyEuler})
    return false
end # due to caching
