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
struct JVODE{bType, aType} <: OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm
    algorithm::Symbol
    bias1::bType
    bias2::bType
    bias3::bType
    addon::aType
end

function JVODE(
        algorithm = :Adams; bias1 = 6, bias2 = 6, bias3 = 10,
        addon = 1 // 10^6
    )
    return JVODE(algorithm, bias1, bias2, bias3, addon)
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
