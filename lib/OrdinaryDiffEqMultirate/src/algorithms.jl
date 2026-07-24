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
    # Allocate dedicated initdt rate-scratch buffers (see Tsit5); default off keeps
    # the cache footprint identical to the historical layout.
    preallocate_initdt_buffers::Bool = false
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
    preallocate_initdt_buffers::Bool = false
end

@doc generic_solver_docstring(
    "Multirate Infinitesimal GARK — explicit midpoint (MRI-GARK-ERK22a).

Solves a SplitODE `du/dt = f1(u,t) + f2(u,t)` where `f1` is the fast component
and `f2` is the slow component (SciML convention). 2-stage, 2nd-order explicit
MRI-GARK method (Sandu 2019); each stage integrates the fast component over a
sub-interval with a slow-coupling, using `m` explicit-midpoint inner micro-steps.
Embedded 1st-order error estimate.",
    "MRIGARKERK22a",
    "Multirate infinitesimal GARK explicit method.",
    """@article{sandu2019class,
    title={A class of multirate infinitesimal {GARK} methods},
    author={Sandu, Adrian},
    journal={SIAM Journal on Numerical Analysis},
    volume={57},
    number={5},
    pages={2300--2327},
    year={2019}}""",
    """
    - `m`: number of inner midpoint micro-steps per stage.
    """,
    """
    m::Int,
    """
)
Base.@kwdef struct MRIGARKERK22a <: OrdinaryDiffEqAdaptiveAlgorithm
    m::Int
    preallocate_initdt_buffers::Bool = false
end

@doc generic_solver_docstring(
    "Multirate Infinitesimal GARK — explicit trapezoidal (MRI-GARK-ERK22b).

Solves a SplitODE `du/dt = f1(u,t) + f2(u,t)` where `f1` is the fast component
and `f2` is the slow component (SciML convention). 2nd-order explicit MRI-GARK
method (the `c₂ = 1` trapezoidal member of the ERK22 family, Sandu 2019), with
a pure slow-correction final stage. Inner fast micro-ODE uses `m` explicit-midpoint
micro-steps. Embedded 1st-order error estimate.",
    "MRIGARKERK22b",
    "Multirate infinitesimal GARK explicit method.",
    """@article{sandu2019class,
    title={A class of multirate infinitesimal {GARK} methods},
    author={Sandu, Adrian},
    journal={SIAM Journal on Numerical Analysis},
    volume={57},
    number={5},
    pages={2300--2327},
    year={2019}}""",
    """
    - `m`: number of inner midpoint micro-steps per stage.
    """,
    """
    m::Int,
    """
)
Base.@kwdef struct MRIGARKERK22b <: OrdinaryDiffEqAdaptiveAlgorithm
    m::Int
    preallocate_initdt_buffers::Bool = false
end

@doc generic_solver_docstring(
    "Multirate Infinitesimal GARK — 3rd-order explicit (MRI-GARK-ERK33a).

Solves a SplitODE `du/dt = f1(u,t) + f2(u,t)` where `f1` is the fast component
and `f2` is the slow component (SciML convention). 3-stage, 3rd-order explicit
MRI-GARK method of Sandu 2019; each stage integrates the fast component over a
sub-interval with a τ-dependent slow-coupling, using `m` explicit-RK3 inner
micro-steps. Embedded 2nd-order error estimate.",
    "MRIGARKERK33a",
    "Multirate infinitesimal GARK explicit method.",
    """@article{sandu2019class,
    title={A class of multirate infinitesimal {GARK} methods},
    author={Sandu, Adrian},
    journal={SIAM Journal on Numerical Analysis},
    volume={57},
    number={5},
    pages={2300--2327},
    year={2019}}""",
    """
    - `m`: number of inner RK3 micro-steps per stage.
    """,
    """
    m::Int,
    """
)
Base.@kwdef struct MRIGARKERK33a <: OrdinaryDiffEqAdaptiveAlgorithm
    m::Int
end

@doc generic_solver_docstring(
    "Multirate Infinitesimal GARK — 4th-order explicit (MRI-GARK-ERK45a).

Solves a SplitODE `du/dt = f1(u,t) + f2(u,t)` where `f1` is the fast component
and `f2` is the slow component (SciML convention). 5-stage, 4th-order explicit
MRI-GARK method of Sandu 2019; each stage integrates the fast component over a
sub-interval with a τ-dependent slow-coupling, using `m` explicit-RK4 inner
micro-steps. Embedded 3rd-order error estimate.",
    "MRIGARKERK45a",
    "Multirate infinitesimal GARK explicit method.",
    """@article{sandu2019class,
    title={A class of multirate infinitesimal {GARK} methods},
    author={Sandu, Adrian},
    journal={SIAM Journal on Numerical Analysis},
    volume={57},
    number={5},
    pages={2300--2327},
    year={2019}}""",
    """
    - `m`: number of inner RK4 micro-steps per stage.
    """,
    """
    m::Int,
    """
)
Base.@kwdef struct MRIGARKERK45a <: OrdinaryDiffEqAdaptiveAlgorithm
    m::Int
end

@doc generic_solver_docstring(
    "Multirate Infinitesimal GARK — 2nd-order solve-decoupled implicit (MRI-GARK-IRK21a).

Solves a SplitODE `du/dt = f1(u,t) + f2(u,t)` where `f1` is the fast component
and `f2` is the slow component (SciML convention). 2nd-order solve-decoupled
implicit MRI-GARK method of Sandu 2019: the fast component is integrated over the
step with `m` explicit-midpoint inner micro-steps and a frozen slow forcing, then
a final `Δc = 0` stage applies an implicit trapezoidal correction that requires a
nonlinear solve in the slow rate `f2`. Suited to problems whose slow component is
itself stiff — the case explicit MRI-GARK cannot handle.",
    "MRIGARKIRK21a",
    "Multirate infinitesimal GARK solve-decoupled implicit method.",
    """@article{sandu2019class,
    title={A class of multirate infinitesimal {GARK} methods},
    author={Sandu, Adrian},
    journal={SIAM Journal on Numerical Analysis},
    volume={57},
    number={5},
    pages={2300--2327},
    year={2019}}""",
    """
    - `m`: number of inner midpoint micro-steps for the fast stage.
    """,
    """
    m::Int,
    """
)
struct MRIGARKIRK21a{AD, F, F2, CJ} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm
    m::Int
    linsolve::F
    nlsolve::F2
    autodiff::AD
    concrete_jac::CJ
end

function MRIGARKIRK21a(;
        m::Int, autodiff = AutoForwardDiff(), concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton()
    )
    autodiff = _fixup_ad(autodiff)
    return MRIGARKIRK21a(m, linsolve, nlsolve, autodiff, _unwrap_val(concrete_jac))
end

@doc generic_solver_docstring(
    "Multirate Infinitesimal GARK — 3rd-order solve-decoupled implicit (MRI-GARK-ESDIRK34a).

Solves a SplitODE `du/dt = f1(u,t) + f2(u,t)` where `f1` is the fast component
and `f2` is the slow component (SciML convention). 6-stage, 3rd-order
solve-decoupled implicit MRI-GARK method of Sandu 2019: the fast component is
integrated over three sub-intervals with `m` explicit-RK3 inner micro-steps, each
followed by a `Δc = 0` stage that applies an ESDIRK slow correction requiring a
nonlinear solve in the slow rate `f2`. The three implicit stages share the common
diagonal coefficient, so a single Jacobian/`W` is reused. Suited to problems whose
slow component is itself stiff — the case explicit MRI-GARK cannot handle.",
    "MRIGARKESDIRK34a",
    "Multirate infinitesimal GARK solve-decoupled implicit method.",
    """@article{sandu2019class,
    title={A class of multirate infinitesimal {GARK} methods},
    author={Sandu, Adrian},
    journal={SIAM Journal on Numerical Analysis},
    volume={57},
    number={5},
    pages={2300--2327},
    year={2019}}""",
    """
    - `m`: number of inner RK3 micro-steps per fast stage.
    """,
    """
    m::Int,
    """
)
struct MRIGARKESDIRK34a{AD, F, F2, CJ} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm
    m::Int
    linsolve::F
    nlsolve::F2
    autodiff::AD
    concrete_jac::CJ
end

function MRIGARKESDIRK34a(;
        m::Int, autodiff = AutoForwardDiff(), concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton()
    )
    autodiff = _fixup_ad(autodiff)
    return MRIGARKESDIRK34a(m, linsolve, nlsolve, autodiff, _unwrap_val(concrete_jac))
end

@doc generic_solver_docstring(
    "Multirate Infinitesimal Step (MIS).

Solves a split ODE of the form `du/dt = f1(u,t) + f2(u,t)` where `f1` is the
fast component and `f2` is the slow component (SciML convention). One MIS
step advances through the tableau's outer stages; each stage solves a small
modified fast ODE on `τ ∈ [0, d_i · dt]` whose right-hand side is the original
fast rate `f1` plus a constant offset built from prior stages' slow tendencies
and α/γ corrections. The inner ODE is approximated with `m · d_i` explicit
midpoint (RK2) micro-steps. Embedded error estimate is `Y_s − Y_{s-1}`.

Currently provides the 2nd-order, 4-stage MIS2(4,2) tableau of
Wensch–Knoth–Galant (BIT 2009).",
    "MIS",
    "Multirate infinitesimal step method.",
    """@article{wensch2009multirate,
    title={Multirate infinitesimal step methods for atmospheric flow simulation},
    author={Wensch, J{\\\"o}rg and Knoth, Oswald and Galant, A},
    journal={BIT Numerical Mathematics},
    volume={49},
    number={2},
    pages={449--473},
    year={2009}}""",
    """
    - `m`: nominal number of inner midpoint micro-steps per `dt`; each stage takes
        `max(1, ceil(m · d_i))` micro-steps.
    """,
    """
    m::Int,
    """
)
Base.@kwdef struct MIS <: OrdinaryDiffEqAdaptiveAlgorithm
    m::Int
end
