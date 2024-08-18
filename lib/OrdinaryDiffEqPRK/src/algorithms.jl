@doc generic_solver_docstring("A 5 parallel, 2 processor method of 5th order.
        ",
        "KuttaPRK2p5",
        "Explicit Runge-Kutta Method",
        "REFERENCE TBD",
        "- `threading`: TBD",
        "threading = true,")
Base.@kwdef struct KuttaPRK2p5{TO} <: OrdinaryDiffEqAlgorithm
    threading::TO = true
end
