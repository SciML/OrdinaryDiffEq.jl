# RK methods

struct ExplicitRK{TabType} <: OrdinaryDiffEqAdaptiveAlgorithm
    tableau::TabType
end
ExplicitRK(; tableau = ODE_DEFAULT_TABLEAU) = ExplicitRK(tableau)

TruncatedStacktraces.@truncate_stacktrace ExplicitRK