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
alg_order(alg::ESDIRK325L2SA) = 3
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

function isesdirk(
        alg::Union{
            KenCarp3, KenCarp4, KenCarp47, KenCarp5, KenCarp58,
            Kvaerno3, Kvaerno4, Kvaerno5, ESDIRK325L2SA, ESDIRK437L2SA,
            ESDIRK54I8L2SA, ESDIRK436L2SA2, ESDIRK547L2SA2,
            ESDIRK659L2SA, CFNLIRK3,
        }
    )
    return true
end

# Methods with an explicit first stage recover `z₁` from the FSAL derivative
# rather than solving for it, which under a mass matrix needs `M⁻¹`. That is
# done elementwise (see `_z_from_dtf`), so these are restricted to diagonal `M`;
# Core's `solve` throws a pointed error for anything else. The implicit-first-stage
# SDIRKs solve for `z₁` and keep general mass-matrix support, and so does Trapezoid,
# whose own step folds `M * uprev` into the stage equation instead of recovering `z₁`.
function only_diagonal_mass_matrix(
        ::Union{
            KenCarp3, KenCarp4, KenCarp47, KenCarp5, KenCarp58,
            Kvaerno3, Kvaerno4, Kvaerno5,
            ESDIRK325L2SA, ESDIRK436L2SA2, ESDIRK437L2SA, ESDIRK547L2SA2,
            ESDIRK54I8L2SA, ESDIRK659L2SA,
            TRBDF2, CFNLIRK3,
            ARS222, ARS232, ARS343, ARS443, BHR553,
        }
    )
    return true
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
issplit(alg::ARS343) = true
alg_order(alg::ARS343) = 3
isesdirk(alg::ARS343) = true
issplit(alg::ARS222) = true
alg_order(alg::ARS222) = 2
isesdirk(alg::ARS222) = true
issplit(alg::ARS232) = true
alg_order(alg::ARS232) = 2
isesdirk(alg::ARS232) = true
issplit(alg::ARS443) = true
alg_order(alg::ARS443) = 3
isesdirk(alg::ARS443) = true
issplit(alg::IMEXSSP222) = true
alg_order(alg::IMEXSSP222) = 2
isesdirk(alg::IMEXSSP222) = true
issplit(alg::IMEXSSP2322) = true
alg_order(alg::IMEXSSP2322) = 2
isesdirk(alg::IMEXSSP2322) = true
issplit(alg::IMEXSSP3332) = true
alg_order(alg::IMEXSSP3332) = 2
isesdirk(alg::IMEXSSP3332) = true
issplit(alg::IMEXSSP3433) = true
alg_order(alg::IMEXSSP3433) = 3
isesdirk(alg::IMEXSSP3433) = true
issplit(alg::BHR553) = true
alg_order(alg::BHR553) = 3
isesdirk(alg::BHR553) = true

# Per-stage Newton-seed strategy. Every SDIRK/ESDIRK algorithm in this module
# carries a `predictor::Predictor.T` field, so the `alg.predictor` access is
# usually direct. The `hasproperty` fallback is for downstream algorithms that
# reuse `ESDIRKIMEXCache` (e.g. OrdinaryDiffEqBDF's `ABDF2`, which uses the
# Implicit Euler tableau as a starter step via `cache.eulercache`) without
# carrying a `predictor` field of their own. Those callers report `Trivial`, but
# the stage-1 seed in `generic_imex_perform_step.jl` re-checks `hasproperty` to
# give them the linear bootstrap seed that preserves their convergence order,
# rather than the zero seed a genuine `Trivial` request selects.
_predictor(alg) = hasproperty(alg, :predictor) ? alg.predictor : Predictor.Trivial

# The interpolant predictors use the Hermite power form; a method with a custom
# interpolant should override this to false to fall back to the full extrapolant.
_uses_hermite_interp(alg) = true

# Deprecated `extrapolant` Symbol -> `Predictor` enum mapping.
function _resolve_predictor(predictor::Predictor.T, extrapolant)
    extrapolant === nothing && return predictor
    Base.depwarn(
        "The `extrapolant` keyword is deprecated; use `predictor` (a `Predictor` enum value).",
        :extrapolant
    )
    extrapolant isa Predictor.T && return extrapolant
    extrapolant === :constant && return Predictor.Trivial
    extrapolant === :linear && return Predictor.Linear
    extrapolant === :interpolant && return Predictor.MaxOrder
    return predictor
end
