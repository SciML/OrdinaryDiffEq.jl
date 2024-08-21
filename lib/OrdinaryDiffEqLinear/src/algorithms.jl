# Linear Methods

for Alg in [
    :MagnusMidpoint,
    :MagnusLeapfrog,
    :LieEuler,
    :MagnusGauss4,
    :MagnusNC6,
    :MagnusGL6,
    :MagnusGL8,
    :MagnusNC8,
    :MagnusGL4,
    :RKMK2,
    :RKMK4,
    :LieRK4,
    :CG2,
    :CG3,
    :CG4a
]
    @eval begin
        @doc generic_solver_docstring("description TBD",
            $(string(Alg)),
            "Semilinear ODE solver",
            "ref TBD",
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

@doc generic_solver_docstring("description TBD",
    "MagnusAdapt4",
    "Semilinear ODE solver",
    "ref TBD", "", "")
struct MagnusAdapt4 <: OrdinaryDiffEqAdaptiveAlgorithm end

@doc generic_solver_docstring("description TBD",
    "CayleyEuler",
    "Semilinear ODE solver",
    "ref TBD", "", "")
struct CayleyEuler <: OrdinaryDiffEqAlgorithm end

@doc generic_solver_docstring("description TBD",
    "LinearExponential",
    "Semilinear ODE solver",
    "ref TBD",
    """
    - `krylov`: TBD
    - `m`: TBD
    - `iop`: TBD
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
