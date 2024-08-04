for Alg in [:LawsonEuler, :NorsettEuler, :ETDRK2, :ETDRK3, :ETDRK4, :HochOst4]
    """
    Hochbruck, Marlis, and Alexander Ostermann. “Exponential Integrators.” Acta
      Numerica 19 (2010): 209–286. doi:10.1017/S0962492910000048.
    """
    @eval struct $Alg{CS, AD, FDT, ST, CJ} <:
                 OrdinaryDiffEqExponentialAlgorithm{CS, AD, FDT, ST, CJ}
        krylov::Bool
        m::Int
        iop::Int
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

for Alg in [:Exprb32, :Exprb43]
    @eval struct $Alg{CS, AD, FDT, ST, CJ} <:
                 OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD, FDT, ST, CJ}
        m::Int
        iop::Int
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
for Alg in [:Exp4, :EPIRK4s3A, :EPIRK4s3B, :EPIRK5s3, :EXPRB53s3, :EPIRK5P1, :EPIRK5P2]
    @eval struct $Alg{CS, AD, FDT, ST, CJ} <:
                 OrdinaryDiffEqExponentialAlgorithm{CS, AD, FDT, ST, CJ}
        adaptive_krylov::Bool
        m::Int
        iop::Int
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