# IMEX Multistep methods

@doc generic_solver_docstring("Crank-Nicholson Adams-Bashforth 2.",
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
    publisher={Wiley Online Library}}", "", "")
struct CNAB2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    autodiff::AD
end

function CNAB2(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    CNAB2{
        _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(
        linsolve,
        nlsolve,
        precs,
        extrapolant,
        AD_choice)
end

@doc generic_solver_docstring("Crank-Nicholson Leapfrong 2.",
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
    publisher={Elsevier}}", "", "")
struct CNLF2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    autodiff::AD
end
function CNLF2(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    CNLF2{
        _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(
        linsolve,
        nlsolve,
        precs,
        extrapolant,
        AD_choice)
end
