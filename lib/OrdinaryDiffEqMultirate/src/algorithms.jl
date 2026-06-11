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

@doc generic_solver_docstring(
    "Multirate Adams–Bashforth (MRAB).

Solves a split ODE of the form `du/dt = f1(u,t) + f2(u,t)` where `f1` is the
fast component and `f2` is the slow component (SciML convention). The slow rate
`f2` is evaluated once per macro interval and held constant across the `m` fast
substeps; the fast rate `f1` is advanced with a `k`-step explicit Adams–Bashforth
formula on each fast substep, using the combined rate `f1 + f2_frozen` for the
linear-multistep update. The first `k - 1` fast substeps of the first macro
interval bootstrap the history with explicit Euler.

Order `k` for the fast component; globally first-order limited by the frozen
slow rate (same starting order as MREEF before Richardson extrapolation).",
    "MRAB",
    "Multirate explicit linear-multistep method.",
    """@techreport{sandu2007multirate,
    title={Multirate explicit Adams methods for time integration of conservation laws},
    author={Sandu, Adrian and Constantinescu, Emil M},
    year={2007},
    institution={Virginia Polytechnic Institute and State University},
    number={CS-TR-07-30}}""",
    """
    - `k`: Adams–Bashforth order, `1 ≤ k ≤ 5`. Default is `2`.
    - `m`: number of fast substeps per macro interval. Default is `4`.
    """,
    """
    k::Int = 2,
    m::Int = 4,
    """
)
Base.@kwdef struct MRAB <: OrdinaryDiffEqAdaptiveAlgorithm
    k::Int = 2
    m::Int = 4
end
