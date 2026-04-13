@doc generic_solver_docstring(
    "Multirate Richardson Extrapolation with Euler as the base method (MREEF).

Solves a split ODE of the form `du/dt = f1(u,t) + f2(u,t)` where `f1` is the
fast component and `f2` is the slow component (SciML convention). The slow rate
`f2` is frozen over each macro interval and the fast rate `f1` is integrated
with `m` explicit Euler substeps. Aitken–Neville Richardson extrapolation over
`order` base solutions is then applied to boost accuracy.",
    "MREEF",
    "Multirate explicit method.",
    """@article{engstrom2009multirate,
    title={Multirate explicit Adams methods for time integration of conservation laws},
    author={Engstr{\\\"o}m, C and Ferm, L and L{\\\"o}tstedt, P and Sj{\\\"o}green, B},
    year={2009}}""",
    """
    - `m`: number of fast substeps per macro interval. Default is `4`.
    - `order`: extrapolation order (number of base solutions). Default is `4`.
    - `seq`: subdivision sequence, `:harmonic` (default) or `:romberg`.
    """,
    """
    m::Int = 4,
    order::Int = 4,
    seq::Symbol = :harmonic,
    """
)
Base.@kwdef struct MREEF <: OrdinaryDiffEqAdaptiveAlgorithm
    m::Int = 4
    order::Int = 4
    seq::Symbol = :harmonic
end
