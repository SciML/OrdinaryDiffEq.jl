# Linear Methods
REF1 = """
@article{celledoni2014introduction,
  title={An introduction to Lie group integrators--basics, new developments and applications},
  author={Celledoni, Elena and Marthinsen, H{\aa}kon and Owren, Brynjulf},
  journal={Journal of Computational Physics},
  volume={257},
  pages={1040--1061},
  year={2014},
  publisher={Elsevier}
}
"""
REF2 = """
@article{crouch1993numerical,
  title={Numerical integration of ordinary differential equations on manifolds},
  author={Crouch, Peter E and Grossman, R},
  journal={Journal of Nonlinear Science},
  volume={3},
  pages={1--33},
  year={1993},
  publisher={Springer}
}
"""
REF3 = """
@article{blanes2000improved,
  title={Improved high order integrators based on the Magnus expansion},
  author={Blanes, Sergio and Casas, Fernando and Ros, Javier},
  journal={BIT Numerical Mathematics},
  volume={40},
  number={3},
  pages={434--450},
  year={2000},
  publisher={Springer}
}
"""
REF4 = """
@article{blanes2009magnus,
  title={The Magnus expansion and some of its applications},
  author={Blanes, Sergio and Casas, Fernando and Oteo, Jose-Angel and Ros, Jos{\'e}},
  journal={Physics reports},
  volume={470},
  number={5-6},
  pages={151--238},
  year={2009},
  publisher={Elsevier}
}
"""
REF5 = """
@article{hairer2011solving,
  title={Solving differential equations on manifolds},
  author={Hairer, Ernst},
  journal={Lecture notes},
  year={2011}
}
"""
REF6 = """
@article{jackiewicz2000construction,
  title={Construction of Runge--Kutta methods of Crouch--Grossman type of high order},
  author={Jackiewicz, Zdzislaw and Marthinsen, Arne and Owren, Brynjulf},
  journal={Advances in Computational Mathematics},
  volume={13},
  pages={405--415},
  year={2000},
  publisher={Springer}
}
"""

for (Alg, Description, Ref) in [
        (
            :MagnusMidpoint, "Second order Magnus Midpoint method.",
            "https://joshuagoings.com/2017/06/15/magnus/",
        ),
        (
            :MagnusLeapfrog, "Second order Magnus Leapfrog method.",
            "https://joshuagoings.com/2017/06/15/magnus/",
        ),
        (:LieEuler, "description", REF1),
        (
            :MagnusGauss4,
            "Fourth order Magnus method approximated using a two stage Gauss quadrature.",
            REF5,
        ),
        (
            :MagnusNC6,
            "Sixth order Magnus method approximated using Newton-Cotes quadrature.", REF3,
        ),
        (
            :MagnusGL6,
            "Sixth order Magnus method approximated using Gauss-Legendre quadrature.", REF3,
        ),
        (
            :MagnusGL8,
            "Eighth order Magnus method approximated using Newton-Cotes quadrature.", REF3,
        ),
        (
            :MagnusNC8,
            "Eighth order Magnus method approximated using Gauss-Legendre quadrature.", REF3,
        ),
        (
            :MagnusGL4,
            "Fourth order Magnus method approximated using Gauss-Legendre quadrature.", REF4,
        ),
        (:RKMK2, "Second order Runge–Kutta–Munthe-Kaas method.", REF1),
        (:RKMK4, "Fourth order Runge–Kutta–Munthe-Kaas method.", REF1),
        (:LieRK4, "Fourth order Lie Runge-Kutta method.", REF1),
        (:CG2, "Second order Crouch–Grossman method.", REF1),
        (:CG3, "Third order Crouch-Grossman method.", REF2),
        (:CG4a, " Fourth order Crouch-Grossman method.", REF6),
    ]
    @eval begin
        @doc generic_solver_docstring(
            $Description,
            $(string(Alg)),
            "Semilinear ODE solver",
            $Ref,
            """
            - `krylov`: Determines whether Krylov approximation or operator caching is used, the latter only available for semilinear problems.
                `krylov=true` is much faster for larger systems and is thus recommended whenever there are >100 ODEs.
            - `m`: Controls the size of Krylov subspace.
            - `iop`: If not zero, determines the length of the incomplete orthogonalization procedure (IOP).
                    Note that if the linear operator/Jacobian is hermitian, then the Lanczos algorithm will always be used and the IOP setting is ignored.
            """,
            """
            krylov = false,
            m = 30,
            iop = 0,
            """
        )
        struct $Alg <: OrdinaryDiffEqLinearExponentialAlgorithm
            krylov::Bool
            m::Int
            iop::Int
        end
    end
    @eval $Alg(; krylov = false, m = 30, iop = 0) = $Alg(krylov, m, iop)
end

@doc generic_solver_docstring(
    "Fourth Order Adaptive Magnus method.",
    "MagnusAdapt4",
    "Semilinear ODE solver",
    "@article{li2008adaptive,
    title={Adaptive explicit Magnus numerical method for nonlinear dynamical systems},
    author={Li, Wen-cheng and Deng, Zi-chen},
    journal={Applied Mathematics and Mechanics},
    volume={29},
    number={9},
    pages={1111--1118},
    year={2008},
    publisher={Springer}}", "", ""
)
struct MagnusAdapt4 <: OrdinaryDiffEqAdaptiveAlgorithm end

@doc generic_solver_docstring(
    "First order method using Cayley transformations.",
    "CayleyEuler",
    "Semilinear ODE solver",
    "@article{iserles2000lie,
    title={Lie-group methods},
    author={Iserles, Arieh and Munthe-Kaas, Hans Z and Norsett, Syvert P and Zanna, Antonella},
    journal={Acta numerica},
    volume={9},
    pages={215--365},
    year={2000},
    publisher={Cambridge University Press}}", "", ""
)
struct CayleyEuler <: OrdinaryDiffEqAlgorithm end

@doc generic_solver_docstring(
    "Exact solution formula for linear, time-independent problems.",
    "LinearExponential",
    "Semilinear ODE solver",
    "@book{strogatz2018nonlinear,
    title={Nonlinear dynamics and chaos: with applications to physics, biology, chemistry, and engineering},
    author={Strogatz, Steven H},
    year={2018},
    publisher={CRC press}}",
    """
    - `krylov`:
        - `:off`: cache the operator beforehand. Requires Matrix(A) method defined for the operator A.
        - `:simple`: uses simple Krylov approximations with fixed subspace size m.
        - `:adaptive`: uses adaptive Krylov approximations with internal timestepping.
    - `m`: Controls the size of Krylov subspace if `krylov=:simple`, and the initial subspace size if `krylov=:adaptive`.
    - `iop`: If not zero, determines the length of the incomplete orthogonalization procedure
            Note that if the linear operator/Jacobian is hermitian,
            then the Lanczos algorithm will always be used and the IOP setting is ignored.
    """,
    """
    krylov = :off,
    m = 10,
    iop = 0,
    """
)
struct LinearExponential <:
    OrdinaryDiffEqExponentialAlgorithm{1, false, Val{:forward}, Val{true}, nothing}
    krylov::Symbol
    m::Int
    iop::Int
end
LinearExponential(; krylov = :off, m = 10, iop = 0) = LinearExponential(krylov, m, iop)
