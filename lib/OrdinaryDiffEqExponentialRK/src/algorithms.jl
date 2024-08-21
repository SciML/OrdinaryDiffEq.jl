REF1 = """
Hochbruck, Marlis, and Alexander Ostermann. “Exponential Integrators.” Acta
  Numerica 19 (2010): 209–286. doi:10.1017/S0962492910000048.
"""
for (Alg, Description, Ref) in [
            (:LawsonEuler, "First order exponential Euler scheme.", REF1),
            (:NorsettEuler, "First order exponential-RK scheme. Alias: `ETD1`", REF1),
            (:ETDRK2, "2nd order exponential-RK scheme.", REF1),
            (:ETDRK3, "3rd order exponential-RK scheme.", REF1),
            (:ETDRK4, "4th order exponential-RK scheme", REF1),
            (:HochOst4, "4th order exponential-RK scheme with stiff order 4.", REF1),
            ]

    @eval begin
            @doc generic_solver_docstring($Description,
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
                """)
            struct $Alg{CS, AD, FDT, ST, CJ} <:
                 OrdinaryDiffEqExponentialAlgorithm{CS, AD, FDT, ST, CJ}
                    krylov::Bool
                    m::Int
                    iop::Int
            end
        end
    @eval function $Alg(; krylov = false, m = 30, iop = 0, autodiff = true,
            standardtag = Val{true}(), concrete_jac = nothing,
            chunk_size = Val{0}(),
            diff_type = Val{:forward})
        $Alg{_unwrap_val(chunk_size), _unwrap_val(autodiff),
            diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(krylov,
            m,
            iop)
    end
end

const ETD1 = NorsettEuler # alias

for (Alg, Description, Ref) in [
            (:Exprb32, "3rd order adaptive Exponential-Rosenbrock scheme.", "ref TBD"),
            (:Exprb43, "4th order adaptive Exponential-Rosenbrock scheme.", "ref TBD")]
    @eval begin @doc generic_solver_docstring($Description,
                        $(string(Alg)),
                        "Semilinear ODE solver",
                        $Ref,
                        """
                        - `m`: Controls the size of Krylov subspace.
                        - `iop`: If not zero, determines the length of the incomplete orthogonalization procedure (IOP).
                            Note that if the linear operator/Jacobian is hermitian, then the Lanczos algorithm will always be used and the IOP setting is ignored.
                        """,
                        """
                        m = 30,
                        iop = 0,
                        """)
                struct $Alg{CS, AD, FDT, ST, CJ} <:
                     OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD, FDT, ST, CJ}
                    m::Int
                    iop::Int
                end
    end
    @eval function $Alg(; m = 30, iop = 0, autodiff = true, standardtag = Val{true}(),
            concrete_jac = nothing, chunk_size = Val{0}(),
            diff_type = Val{:forward})
        $Alg{_unwrap_val(chunk_size), _unwrap_val(autodiff),
            diff_type, _unwrap_val(standardtag),
            _unwrap_val(concrete_jac)}(m,
            iop)
    end
end
for (Alg, Description, Ref) in [
    (:Exp4, "4th order EPIRK scheme.", "REF TBD")
    (:EPIRK4s3A, "4th order EPIRK scheme with stiff order 4.", "REF TBD")
    (:EPIRK4s3B, "4th order EPIRK scheme with stiff order 4.", "REF TBD")
    (:EPIRK5s3, "5th order “horizontal” EPIRK scheme with stiff order 5. Broken.", "REF TBD")
    (:EXPRB53s3, "5th order EPIRK scheme with stiff order 5.", "REF TBD")
    (:EPIRK5P1, "5th order EPIRK scheme", "REF TBD")
    (:EPIRK5P2, "5th order EPIRK scheme", "REF TBD")
    ]
    @eval begin
        @doc generic_solver_docstring($Description,
                $(string(Alg)),
                "Semilinear ODE solver",
                $Ref,
                """
                - `adaptive_krylov`: Determines if the adaptive Krylov algorithm with timestepping of Neisen & Wright is used.
                - `m`: Controls the size of Krylov subspace.
                        - `iop`: If not zero, determines the length of the incomplete orthogonalization procedure (IOP).
                            Note that if the linear operator/Jacobian is hermitian, then the Lanczos algorithm will always be used and the IOP setting is ignored.
                """,
                """
                adaptive_krylov = true,
                m = 30,
                iop = 0,
                """)
            struct $Alg{CS, AD, FDT, ST, CJ} <:
                        OrdinaryDiffEqExponentialAlgorithm{CS, AD, FDT, ST, CJ}
                adaptive_krylov::Bool
                m::Int
                iop::Int
            end
    end
    @eval function $Alg(; adaptive_krylov = true, m = 30, iop = 0, autodiff = true,
            standardtag = Val{true}(), concrete_jac = nothing,
            chunk_size = Val{0}(), diff_type = Val{:forward})
        $Alg{_unwrap_val(chunk_size), _unwrap_val(autodiff), diff_type,
            _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(adaptive_krylov,
            m,
            iop)
    end
end

"""
ETD2: Exponential Runge-Kutta Method
Second order Exponential Time Differencing method (in development).
"""
struct ETD2 <:
       OrdinaryDiffEqExponentialAlgorithm{0, false, Val{:forward}, Val{true}, nothing} end
