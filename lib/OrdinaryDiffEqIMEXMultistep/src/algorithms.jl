# IMEX Multistep methods

@doc generic_solver_docstring(
    "Crank-Nicolson Adams Bashforth Order 2 (fixed time step)",
    "CNAB2",
    "IMEX Multistep method.",
    "@article{jorgenson2014unconditional,
    title={Unconditional stability of a Crank-Nicolson Adams-Bashforth 2 numerical method},
    author={JORGENSON, ANDREW D},
    journal={A (A- C)},
    volume={1},
    number={2},
    pages={1},
    year={2014}}
    @article{he2010numerical,
    title={Numerical implementation of the Crank--Nicolson/Adams--Bashforth scheme for the time-dependent Navier--Stokes equations},
    author={He, Yinnian and Li, Jian},
    journal={International journal for numerical methods in fluids},
    volume={62},
    number={6},
    pages={647--659},
    year={2010},
    publisher={Wiley Online Library}}", "", ""
)
struct CNAB2{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end

function CNAB2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear
    )
    autodiff = _fixup_ad(autodiff)

    return CNAB2(
        linsolve,
        nlsolve,
        extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc generic_solver_docstring(
    "Crank-Nicholson Leapfrong 2.",
    "CNLF2",
    "IMEX Multistep method.",
    "@article{han2020second,
    title={A second order, linear, unconditionally stable, Crank--Nicolson--Leapfrog scheme for phase field models of two-phase incompressible flows},
    author={Han, Daozhi and Jiang, Nan},
    journal={Applied Mathematics Letters},
    volume={108},
    pages={106521},
    year={2020},
    publisher={Elsevier}}
    @article{jiang2015crank,
    title={A Crank--Nicolson Leapfrog stabilization: Unconditional stability and two applications},
    author={Jiang, Nan and Kubacki, Michaela and Layton, William and Moraiti, Marina and Tran, Hoang},
    journal={Journal of Computational and Applied Mathematics},
    volume={281},
    pages={263--276},
    year={2015},
    publisher={Elsevier}}", "", ""
)
struct CNLF2{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end
function CNLF2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear
    )
    autodiff = _fixup_ad(autodiff)

    return CNLF2(
        linsolve,
        nlsolve,
        extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end
