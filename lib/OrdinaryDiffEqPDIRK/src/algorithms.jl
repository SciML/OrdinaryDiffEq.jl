@doc differentiation_rk_docstring(
    "A 2 processor 4th order diagonally non-adaptive implicit method.",
    "PDIRK44",
    "Parallel Diagonally Implicit Runge-Kutta Method.";
    references = """"@article{iserles1990theory,
    title={On the theory of parallel Runge—Kutta methods},
    author={Iserles, Arieh and Norrsett, SP},
    journal={IMA Journal of numerical Analysis},
    volume={10},
    number={4},
    pages={463--488},
    year={1990},
    publisher={Oxford University Press}}""",
    extra_keyword_description = """
    - `nlsolve`: TBD,
    - `extrapolant`: TBD,
    - `thread`: determines whether internal broadcasting on appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when Julia is started with multiple threads.
    """,
    extra_keyword_default = """
    nlsolve = NLNewton(),
    extrapolant = :constant,
    thread = OrdinaryDiffEq.True(),
    """
)
struct PDIRK44{AD, F, F2, P, TO} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    threading::TO
    autodiff::AD
    concrete_jac::Union{Nothing, Bool}
end
function PDIRK44(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant, threading = true
    )
    autodiff = _fixup_ad(autodiff)

    return PDIRK44(
        linsolve, nlsolve, precs,
        extrapolant, threading, autodiff,
        _unwrap_val(concrete_jac)

    )
end
