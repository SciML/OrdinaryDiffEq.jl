alg_extrapolates(alg::ImplicitEuler) = true
alg_extrapolates(alg::Trapezoid) = true
alg_extrapolates(alg::SDIRK22) = true

alg_order(alg::Trapezoid) = 2
alg_order(alg::ImplicitEuler) = 1
alg_order(alg::ImplicitMidpoint) = 2
alg_order(alg::TRBDF2) = 2
alg_order(alg::SSPSDIRK2) = 2
alg_order(alg::SDIRK2) = 2
alg_order(alg::SDIRK22) = 2
alg_order(alg::Kvaerno3) = 3
alg_order(alg::Kvaerno4) = 4
alg_order(alg::Kvaerno5) = 5
alg_order(alg::ESDIRK54I8L2SA) = 5
alg_order(alg::ESDIRK436L2SA2) = 4
alg_order(alg::ESDIRK437L2SA) = 4
alg_order(alg::ESDIRK547L2SA2) = 5
alg_order(alg::ESDIRK659L2SA) = 6
alg_order(alg::KenCarp3) = 3
alg_order(alg::CFNLIRK3) = 3
alg_order(alg::KenCarp4) = 4
alg_order(alg::KenCarp47) = 4
alg_order(alg::KenCarp5) = 5
alg_order(alg::KenCarp58) = 5
alg_order(alg::Cash4) = 4
alg_order(alg::SFSDIRK4) = 4
alg_order(alg::SFSDIRK5) = 4
alg_order(alg::SFSDIRK6) = 4
alg_order(alg::SFSDIRK7) = 4
alg_order(alg::SFSDIRK8) = 4
alg_order(alg::Hairer4) = 4
alg_order(alg::Hairer42) = 4

function isesdirk(alg::Union{KenCarp3, KenCarp4, KenCarp5, KenCarp58,
        Kvaerno3, Kvaerno4, Kvaerno5, ESDIRK437L2SA,
        ESDIRK54I8L2SA, ESDIRK436L2SA2, ESDIRK547L2SA2,
        ESDIRK659L2SA, CFNLIRK3})
    true
end

alg_adaptive_order(alg::Trapezoid) = 1
alg_adaptive_order(alg::ImplicitMidpoint) = 1
alg_adaptive_order(alg::ImplicitEuler) = 0

ssp_coefficient(alg::SSPSDIRK2) = 4

isesdirk(alg::TRBDF2) = true

issplit(alg::KenCarp3) = true
issplit(alg::KenCarp4) = true
issplit(alg::KenCarp47) = true
issplit(alg::KenCarp5) = true
issplit(alg::KenCarp58) = true
issplit(alg::CFNLIRK3) = true