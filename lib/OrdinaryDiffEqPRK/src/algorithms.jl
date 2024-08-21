@doc generic_solver_docstring("A 5 parallel, 2 processor method of 5th order.",
    "KuttaPRK2p5",
    "Explicit Runge-Kutta Method",
    """@article{jackson1995potential,
    title={The potential for parallelism in Runge--Kutta methods. Part 1: RK formulas in standard form},
    author={Jackson, Kenneth R and Norsett, Syvert Paul},
    journal={SIAM journal on numerical analysis},
    volume={32},
    number={1},
    pages={49--82},
    year={1995},
    publisher={SIAM}}""",
    "- `threading`: TBD",
    "threading = true,")
Base.@kwdef struct KuttaPRK2p5{TO} <: OrdinaryDiffEqAlgorithm
    threading::TO = true
end
