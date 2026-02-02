# Adams/BDF methods in Nordsieck forms
@doc generic_solver_docstring(
    """An adaptive 5th order fixed-leading coefficient Adams method in Nordsieck form.
    !!! warning "Experimental"
        `AN5` is experimental, the solver `VCABM` is generally preferred.
    """,
    "AN5",
    "Adaptive step size Adams explicit Method",
    "",
    "",
    ""
)
struct AN5 <: OrdinaryDiffEqAdaptiveAlgorithm end
"""
!!! warning "Experimental"

    `JVODE` is experimental, the solver `VCABM` is generally preferred.
"""
struct JVODE{bType, aType, qType} <: OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm
    algorithm::Symbol
    bias1::bType
    bias2::bType
    bias3::bType
    addon::aType
    qmax::qType
    qsteady_min::qType
    qsteady_max::qType
end

function JVODE(
        algorithm = :Adams; bias1 = 6, bias2 = 6, bias3 = 10,
        addon = 1 // 10^6, qmax = float(10), qsteady_min = float(1), qsteady_max = float(3 // 2)
    )
    return JVODE(algorithm, bias1, bias2, bias3, addon, qmax, qsteady_min, qsteady_max)
end
"""
!!! warning "Experimental"

    `JVODE` is experimental, the solver `VCABM` is generally preferred.
"""
JVODE_Adams(; kwargs...) = JVODE(:Adams; kwargs...)
"""
!!! warning "Experimental"

    `JVODE` is experimental, the solver `FBDF` is generally preferred.
"""
JVODE_BDF(; kwargs...) = JVODE(:BDF; kwargs...)
