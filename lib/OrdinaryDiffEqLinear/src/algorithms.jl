# Linear Methods

for (Alg, Description, Ref) in [
    (:MagnusMidpoint, "Second order Magnus Midpoint method.", "ref TBD"),
    (:MagnusLeapfrog, "Second order Magnus Leapfrog method.", "ref TBD"),
    (:LieEuler, "description", "ref TBD"),
    (:MagnusGauss4, "Fourth order Magnus method approximated using a two stage Gauss quadrature.", "ref TBD"),
    (:MagnusNC6, "Sixth order Magnus method approximated using Newton-Cotes quadrature.", "ref TBD"),
    (:MagnusGL6, "Sixth order Magnus method approximated using Gauss-Legendre quadrature.", "ref TBD"),
    (:MagnusGL8, "Eighth order Magnus method approximated using Newton-Cotes quadrature.", "ref TBD"),
    (:MagnusNC8, "Eighth order Magnus method approximated using Gauss-Legendre quadrature.", "ref TBD"),
    (:MagnusGL4, "Fourth order Magnus method approximated using Gauss-Legendre quadrature.", "ref TBD"),
    (:RKMK2, "Second order Runge–Kutta–Munthe-Kaas method.", "ref TBD"),
    (:RKMK4, "Fourth order Runge–Kutta–Munthe-Kaas method.", "ref TBD"),
    (:LieRK4, "Fourth order Lie Runge-Kutta method.", "ref TBD"),
    (:CG2, "Second order Crouch–Grossman method.", "ref TBD"),
    (:CG3, "Third order Crouch-Grossman method.", "ref TBD"),
    (:CG4a, " Fourth order Crouch-Grossman method.", "ref TBD")]
    @eval begin
        @doc generic_solver_docstring($Description,
            $(string(Alg)),
            "Semilinear ODE solver",
            $Ref,
            """
            - `krylov`: TBD
            - `m`: TBD
            - `iop`: TBD
            """,
            """
            krylov = false,
            m = 30,
            iop = 0,
            """)
        struct $Alg <: OrdinaryDiffEqLinearExponentialAlgorithm
        krylov::Bool
        m::Int
        iop::Int
        end
    end
    @eval $Alg(; krylov = false, m = 30, iop = 0) = $Alg(krylov, m, iop)
end

@doc generic_solver_docstring("Fourth Order Adaptive Magnus method.",
    "MagnusAdapt4",
    "Semilinear ODE solver",
    "ref TBD", "", "")
struct MagnusAdapt4 <: OrdinaryDiffEqAdaptiveAlgorithm end

@doc generic_solver_docstring("First order method using Cayley transformations.",
    "CayleyEuler",
    "Semilinear ODE solver",
    "ref TBD", "", "")
struct CayleyEuler <: OrdinaryDiffEqAlgorithm end

@doc generic_solver_docstring("Exact solution formula for linear, time-independent problems.",
    "LinearExponential",
    "Semilinear ODE solver",
    "ref TBD",
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
    """)
struct LinearExponential <:
       OrdinaryDiffEqExponentialAlgorithm{1, false, Val{:forward}, Val{true}, nothing}
    krylov::Symbol
    m::Int
    iop::Int
end
LinearExponential(; krylov = :off, m = 10, iop = 0) = LinearExponential(krylov, m, iop)
