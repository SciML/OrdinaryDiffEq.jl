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
    """)
struct PDIRK44{CS, AD, F, F2, P, FDT, ST, CJ, TO} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    threading::TO
    autodiff::AD
end
function PDIRK44(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant, threading = true)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    PDIRK44{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(threading)}(linsolve, nlsolve, precs,
        extrapolant, threading, AD_choice)
end
