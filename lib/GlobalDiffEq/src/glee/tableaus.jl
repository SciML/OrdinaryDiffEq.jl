# GLEE (global-error-estimating general linear method) tableaus in the y-ε form
# of Constantinescu (2016), doi:10.1137/15M1014633 (arXiv:1503.05166):
#
#   Y_i    = U[i,1] y⁽ⁿ⁻¹⁾ + U[i,2] ε⁽ⁿ⁻¹⁾ + Δt Σ_j A[i,j] f(Y_j)
#   y⁽ⁿ⁾   = y⁽ⁿ⁻¹⁾ + Δt Σ_j B[1,j] f(Y_j)
#   ε⁽ⁿ⁾   = ε⁽ⁿ⁻¹⁾ + Δt Σ_j B[2,j] f(Y_j)
#
# ε⁽ⁿ⁾ is an asymptotically correct estimate of the global error at t_n, and the
# per-step increment ε⁽ⁿ⁾ − ε⁽ⁿ⁻¹⁾ is an asymptotically correct local error
# estimate for the order-p solution y⁽ⁿ⁾.
struct GLEETableau{MType, M2Type, BType, cType, T}
    A::MType
    U::M2Type
    B::BType
    c::cType
    γ::T
    order::Int
    # Stage 1 equals the propagated solution register (U row [1, 0], zero A row),
    # so its f evaluation can reuse fsalfirst.
    stage1_fsal::Bool
    # Index of a stage that equals the step's end solution y⁽ⁿ⁾ (FSAL-style), so
    # its f evaluation doubles as fsallast; 0 when no such stage exists.
    solution_stage::Int
end

nstages(tab::GLEETableau) = length(tab.c)

function _detect_solution_stage(A, U, B, c)
    for i in 2:length(c)
        isone(U[i, 1]) && iszero(U[i, 2]) || continue
        isone(c[i]) || continue
        all(j -> A[i, j] == B[1, j], 1:(i - 1)) || continue
        all(iszero, @view B[1, i:end]) || continue
        return i
    end
    return 0
end

# Convert a tableau from the paper's y-ỹ form to the y-ε form via Lemma 4.1 with
# T = [1 0; 1 1-γ]: U_yε[:,1] = U_yỹ*1, U_yε[:,2] = (1-γ) U_yỹ[:,2],
# B_yε = T⁻¹ B_yỹ.
function _yytilde_to_yeps(U::AbstractMatrix, B::AbstractMatrix, γ)
    Uyε = hcat(U[:, 1] .+ U[:, 2], (1 - γ) .* U[:, 2])
    Byε = vcat(B[1:1, :], (B[2:2, :] .- B[1:1, :]) ./ (1 - γ))
    return Uyε, Byε
end

function _glee_tableau(::Type{T}, ::Type{T2}, A, U, B, γ, order) where {T, T2}
    c_exact = vec(sum(A, dims = 2))
    stage1_fsal = all(iszero, A[1, :]) && isone(U[1, 1]) && iszero(U[1, 2])
    solution_stage = _detect_solution_stage(A, U, B, c_exact)
    return GLEETableau(
        map(T, A), map(T, U), map(T, B), map(T2, c_exact), T(γ), order,
        stage1_fsal, solution_stage
    )
end

"""
    GLEE23Tableau(T, T2)

Tableau of the 3-stage, order-2 explicit GLEE method of Constantinescu (2016),
eq. (4.6) (PETSc's `TSGLEE23`). `γ = 0`, `B·U` diagonal.
"""
function GLEE23Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    A = [
        0 0 0
        1 0 0
        1 // 4 1 // 4 0
    ]
    U = [
        1 0
        1 10
        1 -1
    ]
    B = [
        1 // 12 1 // 12 5 // 6
        1 // 12 1 // 12 -1 // 6
    ]
    return _glee_tableau(T, T2, A, U, B, 0, 2)
end

"""
    GLEE24Tableau(T, T2)

Tableau of the 4-stage, order-2 explicit GLEE method of Constantinescu (2016),
eq. (A.3) (PETSc's `TSGLEE24`). `γ = 0`; both decoupling conditions (`B·U` and
`B·A·U` diagonal) hold, making it the robust second-order choice for long-time
integration.
"""
function GLEE24Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    A = [
        0 0 0 0
        3 // 4 0 0 0
        1 // 4 29 // 60 0 0
        -21 // 44 145 // 44 -20 // 11 0
    ]
    # y-ỹ form as printed in the paper
    U = [
        0 1
        75 // 58 -17 // 58
        0 1
        0 1
    ]
    B = [
        109 // 275 58 // 75 -37 // 110 1 // 6
        3 // 11 0 75 // 88 -1 // 8
    ]
    Uyε, Byε = _yytilde_to_yeps(U, B, 0 // 1)
    return _glee_tableau(T, T2, A, Uyε, Byε, 0, 2)
end

"""
    GLEE35Tableau(T, T2)

Tableau of the 5-stage, order-3 explicit GLEE method of Constantinescu (2016),
eq. (4.9) (PETSc's `TSGLEE35`). `γ = 0`; both decoupling conditions hold. The
coefficients are exact rationals whose numerators exceed `Int64`, so the
tableau is constructed with `BigInt` rationals before conversion.
"""
function GLEE35Tableau(::Type{T}, ::Type{T2}) where {T, T2}
    R(num, den) = big(num) // big(den)
    Z = R(0, 1)
    A = [
        Z Z Z Z Z;
        R(-2169604947363702313, 24313474998937147335) Z Z Z Z;
        R(46526746497697123895, 94116917485856474137) R(-10297879244026594958, 49199457603717988219) Z Z Z;
        R(23364788935845982499, 87425311444725389446) R(-79205144337496116638, 148994349441340815519) R(40051189859317443782, 36487615018004984309) Z Z;
        R(42089522664062539205, 124911313006412840286) R(-15074384760342762939, 137927286865289746282) R(-62274678522253371016, 125918573676298591413) R(13755475729852471739, 79257927066651693390) Z
    ]
    # y-ỹ form as printed in the paper
    U = [
        R(70820309139834661559, 80863923579509469826) R(10043614439674808267, 80863923579509469826);
        R(161694774978034105510, 106187653640211060371) R(-55507121337823045139, 106187653640211060371);
        R(78486094644566264568, 88171030896733822981) R(9684936252167558413, 88171030896733822981);
        R(65394922146334854435, 84570853840405479554) R(19175931694070625119, 84570853840405479554);
        R(8607282770183754108, 108658046436496925911) R(100050763666313171803, 108658046436496925911)
    ]
    B = [
        R(61546696837458703723, 56982519523786160813) R(-55810892792806293355, 206957624151308356511) R(24061048952676379087, 158739347956038723465) R(3577972206874351339, 7599733370677197135) R(-59449832954780563947, 137360038685338563670);
        R(-9738262186984159168, 99299082461487742983) R(-32797097931948613195, 61521565616362163366) R(42895514606418420631, 71714201188501437336) R(22608567633166065068, 55371917805607957003) R(94655809487476459565, 151517167160302729021)
    ]
    Uyε, Byε = _yytilde_to_yeps(U, B, Z)
    return _glee_tableau(T, T2, A, Uyε, Byε, 0, 3)
end
