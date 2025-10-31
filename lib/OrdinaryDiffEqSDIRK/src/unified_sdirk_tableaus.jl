using StaticArrays

abstract type AbstractTableau{T} end

struct SDIRKTableau{T, T2, S, hasEmbedded, hasAdditiveSplitting, hasExplicit, hasPred} <: AbstractTableau{T}
    A::SMatrix{S, S, T}
    b::SVector{S, T}
    c::SVector{S, T2}
    b_embed::Union{SVector{S, T}, Nothing}
    γ::T
    order::Int
    embedded_order::Int
    is_fsal::Bool
    is_stiffly_accurate::Bool
    is_A_stable::Bool
    is_L_stable::Bool
    predictor_type::Symbol
    A_explicit::Union{SMatrix{S, S, T}, Nothing}
    b_explicit::Union{SVector{S, T}, Nothing}
    c_explicit::Union{SVector{S, T2}, Nothing}
    α_pred::Union{SMatrix{S, S, T2}, Nothing}
    has_spice_error::Bool
end

function SDIRKTableau(A::SMatrix{S, S, T}, b::SVector{S, T}, c::SVector{S, T2}, γ::T,
                      order::Int; b_embed=nothing, embedded_order=0,
                      is_fsal=false, is_stiffly_accurate=false,
                      is_A_stable=true, is_L_stable=false,
                      predictor_type=:default, has_additive_splitting=false,
                      A_explicit=nothing, b_explicit=nothing, c_explicit=nothing,
                      α_pred=nothing, has_spice_error=false) where {S, T, T2}
    
    hasEmbedded = b_embed !== nothing
    hasAdditiveSplitting = has_additive_splitting
    hasExplicit = A_explicit !== nothing
    hasPred = α_pred !== nothing
    SDIRKTableau{T, T2, S, hasEmbedded, hasAdditiveSplitting, hasExplicit, hasPred}(
        A, b, c, b_embed, γ, order, embedded_order,
        is_fsal, is_stiffly_accurate, is_A_stable,
        is_L_stable, predictor_type,
        A_explicit, b_explicit, c_explicit, α_pred, has_spice_error)
end

function TRBDF2Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γ = T2(2 - sqrt(2))
    d = T(1 - sqrt(2) / 2) 
    ω = T(sqrt(2) / 4)
    
    A = @SMatrix [0  0  0;
                  d  d  0;
                  ω  ω  d]
    
    b = @SVector [ω, ω, d]
    c = @SVector [0, γ, 1]
    
    b_embed = @SVector [(1-ω)/3, (3*ω+1)/3, d/3]
    
    α1 = T2(-sqrt(2) / 2)
    α2 = T2(1 + sqrt(2) / 2)
    α_pred = @SMatrix T2[
        0 0 0;
        0 0 0;
        α1 α2 0
    ]
    
    SDIRKTableau(A, b, c, d, 2;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:trbdf2_special,
                 α_pred=α_pred)
end

function ImplicitEulerTableau(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    A = @SMatrix [T(1.0)]
    b = @SVector [T(1.0)]
    c = @SVector [T2(1.0)]
    γ = T(1.0)
    
    SDIRKTableau(A, b, c, γ, 1;
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:trivial,
                 has_spice_error=true)
end

function ImplicitMidpointTableau(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.5)
    γc = T2(0.5)
    A = @SMatrix [γT]
    b = @SVector [T(1.0)]
    c = @SVector [γc]
    
    SDIRKTableau(A, b, c, γT, 2;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:trivial)
end

function TrapezoidTableau(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.5)
    A = @SMatrix [T(0)   T(0);
                  T(0.5) T(0.5)]
    b = @SVector [T(0.5), T(0.5)]
    c = @SVector [T2(0.0), T2(1.0)]
    
    SDIRKTableau(A, b, c, γT, 2;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:default,
                 has_spice_error=true)
end

function SDIRK2Tableau(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(1 - 1/sqrt(2))
    γc = T2(1 - 1/sqrt(2))
    A = @SMatrix [γT        T(0);
                  T(1)-γT   γT]
    b = @SVector [T(1)-γT, γT]
    c = @SVector [γc, T2(1)]
    
    SDIRKTableau(A, b, c, γT, 2;
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:default)
end

function SSPSDIRK2Tableau(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(1 - 1/sqrt(2))
    γc = T2(1 - 1/sqrt(2))
    A = @SMatrix [γT     T(0);
                  T(1)-T(2)*γT   γT]
    b = @SVector [T(0.5), T(0.5)]
    c = @SVector [γc, T2(1)-γc]
    
    SDIRKTableau(A, b, c, γT, 2;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:default)
end

function Kvaerno3Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.4358665215)
    γc = T2(0.4358665215)
    
    A = @SMatrix [γT                       0                        0                      0;
                  T(0.490563388419108)     γT                       0                      0;
                  T(0.073570090080892)     0                        γT                     0;
                  T(0.308809969973036)     T(1.490563388254106)     -T(1.235239879727145)  γT]
    
    b = @SVector [T(0.490563388419108), T(0.073570090080892), T(0.4358665215), T(0.0)]
    c = @SVector [γc, 2γc, T2(1), T2(1)]
    
    b_embed = b - A[4, :]
    
    # Build Hermite-style predictor coefficients for stage 3 from z₁ and z₂
    # θ = c3/c2 over interval [0,c2]
    c2 = 2γc
    θ = c[3] / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γc)
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γc)
    α_pred = @SMatrix T2[
        0 0 0 0;
        0 0 0 0;
        α31 α32 0 0;
        0 0 0 0
    ]

    SDIRKTableau(A, b, c, γT, 3;
                 b_embed=b_embed, embedded_order=2,
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite,
                 α_pred=α_pred)
end

function KenCarp3Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.435866521508459)
    γc = T2(0.435866521508459)
    
    A = @SMatrix [γT                         0                        0                       0;
                  T(0.2576482460664272)      γT                       0                       0;
                  -T(0.09351476757488625)    0                        γT                      0;
                  T(0.18764102434672383)     -T(0.595297473576955)    T(0.9717899277217721)   γT]
    
    b = @SVector [T(2756255671327//12835298489170), 
                  -T(10771552573575//22201958757719),
                  T(9247589265047//10645013368117), 
                  T(2193209047091//5459859503100)]
                  
    c = @SVector [γc, 2γc, T2(0.6), T2(1)]
    
    b_embed = A[4, :]
    
    A_explicit = @SMatrix [0                          0                        0                      0;
                           T(0.871733043016918)      0                        0                      0;
                           T(0.5275890119763004)     T(0.0724109880236996)    0                      0;
                           T(0.3990960076760701)     -T(0.4375576546135194)   T(1.0384616469374492)  0]
    
    b_explicit = @SVector [T(0.18764102434672383), 
                           -T(0.595297473576955),
                           T(0.9717899277217721), 
                           T(0.435866521508459)]
    
    c_explicit = @SVector [T2(0), γc, T2(0.6), T2(1)]
    
    # Build Hermite-style predictor coefficients for stages 3 and 4
    c2 = 2γc
    θ = c[3] / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γc)
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γc)
    θ4 = c[4] / c2 # == 1
    α41 = ((1 + (-4θ4 + 3θ4^2)) + (6θ4 * (1 - θ4) / c2) * γc)
    α42 = ((-2θ4 + 3θ4^2) + (6θ4 * (1 - θ4) / c2) * γc)
    α_pred = @SMatrix T2[
        0 0 0 0;
        0 0 0 0;
        α31 α32 0 0;
        α41 α42 0 0
    ]

    SDIRKTableau(A, b, c, γT, 3;
                 b_embed=b_embed, embedded_order=2,
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:kencarp_additive,
                 has_additive_splitting=true,
                 A_explicit=A_explicit, b_explicit=b_explicit, c_explicit=c_explicit,
                 α_pred=α_pred)
end

function Kvaerno4Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.4358665215)
    γc = T2(0.4358665215)
    
    A = @SMatrix [γT                       0                        0                       0                       0;
                  T(0.490563388419108)     γT                       0                       0                       0;
                  T(0.073570090080892)     0                        γT                      0                       0;
                  T(0.308809969973036)     T(1.490563388254106)     -T(1.235239879727145)   γT                      0;
                  T(0.490563388419108)     T(0.073570090080892)     T(0.4358665215)         T(0.0)                  γT]
    
    b = @SVector [T(0.490563388419108), T(0.073570090080892), T(0.4358665215), T(0.0), T(0.0)]
    c = @SVector [γc, 2γc, T2(1), T2(1), T2(1)]
    
    b_embed = A[5, :]
    
    SDIRKTableau(A, b, c, γT, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function Kvaerno5Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.26)
    γc = T2(0.26)
    
    A = @SMatrix [γT                       0                        0                       0                       0                       0;
                  T(0.13)                  γT                       0                       0                       0                       0;
                  T(0.84079895052208)      T(0.07920104947792)      γT                      0                       0                       0;
                  T(0.619100897516618)     T(0.066593016584582)     T(0.054305985899400)    γT                      0                       0;
                  T(0.258023287184119)     T(0.097741417057132)     T(0.464732297848610)    T(1.179502539939939)    γT                      0;
                  T(0.544974750228521)     T(0.212765981366776)     T(0.164488906111538)    T(0.077770561901165)    T(0.0)                  γT]
    
    b = @SVector [T(0.544974750228521), T(0.212765981366776), T(0.164488906111538), T(0.077770561901165), T(0.0), T(0.0)]
    c = @SVector [γc, T2(0.39), T2(1.0), T2(0.74), T2(1.0), T2(1.0)]
    
    b_embed = A[6, :]
    
    SDIRKTableau(A, b, c, γT, 5;
                 b_embed=b_embed, embedded_order=4,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function KenCarp4Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(1//4)
    γc = T2(1//4)
    
    A = @SMatrix [γT                            0                             0                            0                            0                            0;
                  T(1//2)                       γT                            0                            0                            0                            0;
                  T(83//250)                    T(-13//250)                   γT                           0                            0                            0;
                  T(31//50)                     T(-11//20)                    T(11//20)                    γT                           0                            0;
                  T(17//20)                     T(-1//4)                      T(1//4)                      T(1//2)                      γT                           0;
                  T(755//1728)                  T(755//1728)                  T(-1640//1728)              T(1640//1728)                T(1//4)                      γT]
    
    b = @SVector [T(755//1728), T(755//1728), T(-1640//1728), T(1640//1728), T(1//4), T(0)]
    c = @SVector [γc, T2(3//4), T2(11//20), T2(1//2), T2(1), T2(1)]
    
    b_embed = A[6, :]
    
    A_explicit = @SMatrix [0                             0                             0                            0                            0                            0;
                           T(1//2)                       0                             0                            0                            0                            0;
                           T(13861//62500)               T(6889//62500)                0                            0                            0                            0;
                           T(-116923316275//2393684061468) T(-2731218467317//15368042101831) T(9408046702089//11113171139209) 0                            0                            0;
                           T(-451086348788//2902428689909) T(-2682348792572//7519795681897) T(12662868775082//11960479115383) T(3355817975965//11060851509271) 0                            0;
                           T(647845179188//3216320057751) T(73281519250//8382639484533) T(552539513391//3454668386233) T(3354512671639//8306763924573) T(4040//17871)               0]
    
    b_explicit = @SVector [T(82889//524892), T(0), T(15625//83664), T(69875//102672), T(-2260//8211), T(1//4)]
    c_explicit = @SVector [0, T2(1//2), T2(83//250), T2(31//50), T2(17//20), T2(1)]
    
    SDIRKTableau(A, b, c, γT, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:kencarp_additive,
                 has_additive_splitting=true,
                 A_explicit=A_explicit, b_explicit=b_explicit, c_explicit=c_explicit)
end

function KenCarp47Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.1496590219993)
    γc = T2(0.1496590219993)
    
    A = @SMatrix [γT                            0                             0                            0                            0                            0                            0;
                  T(0.7481896206814)            γT                            0                            0                            0                            0                            0;
                  T(0.2068513093527)            T(0.2931486906473)            γT                           0                            0                            0                            0;
                  T(0.7581896206812)            T(-0.2581896206812)           T(0.25)                      γT                           0                            0                            0;
                  T(0.8765725810946)            T(-0.3765725810946)           T(0.25)                      T(0.25)                      γT                           0                            0;
                  T(1.6274999742127)            T(-1.1274999742127)           T(0.25)                      T(0.25)                      T(0.0)                       γT                           0;
                  T(1.6274999742127)            T(-1.1274999742127)           T(0.25)                      T(0.25)                      T(0.0)                       T(0.0)                       γT]
    
    b = @SVector [T(1.6274999742127), T(-1.1274999742127), T(0.25), T(0.25), T(0.0), T(0.0), T(0.0)]
    c = @SVector [γc, T2(0.8978486300007), T2(0.6496590219993), T2(0.7496590219993), T2(1.1262315616939), T2(1.0), T2(1.0)]
    
    b_embed = A[7, :]
    
    A_explicit = zeros(SMatrix{7,7,T})
    b_explicit = zeros(SVector{7,T})
    c_explicit = c
    
    SDIRKTableau(A, b, c, γT, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:kencarp_additive,
                 has_additive_splitting=true,
                 A_explicit=A_explicit, b_explicit=b_explicit, c_explicit=c_explicit)
end

function KenCarp5Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.2113248654051871)
    γc = T2(0.2113248654051871)
    
    A = @SMatrix [γT                            0                             0                            0                            0                            0                            0                            0;
                  T(0.5)                        γT                            0                            0                            0                            0                            0                            0;
                  T(0.453125)                   T(-0.140625)                  γT                           0                            0                            0                            0                            0;
                  T(0.6828)                     T(-0.2178)                    T(0.3237)                    γT                           0                            0                            0                            0;
                  T(0.6262)                     T(-0.1848)                    T(0.2477)                    T(0.4998)                    γT                           0                            0                            0;
                  T(0.3415)                     T(-0.1219)                    T(0.2502)                    T(0.2502)                    T(0.0686)                    γT                           0                            0;
                  T(0.3415)                     T(-0.1219)                    T(0.2502)                    T(0.2502)                    T(0.0686)                    T(0.0)                       γT                           0;
                  T(0.6262)                     T(-0.1848)                    T(0.2477)                    T(0.4998)                    T(0.2113)                    T(0.0)                       T(0.0)                       γT]
    
    b = @SVector [T(0.3415), T(-0.1219), T(0.2502), T(0.2502), T(0.0686), T(0.0), T(0.0), T(0.0)]
    c = @SVector [γc, T2(0.7113248654051871), T2(0.5226873345948129), T2(0.7887), T2(0.9773126654051871), T2(1.0), T2(1.0), T2(1.0)]
    
    b_embed = A[8, :]
    
    A_explicit = zeros(SMatrix{8,8,T})
    b_explicit = zeros(SVector{8,T})
    c_explicit = c
    
    SDIRKTableau(A, b, c, γT, 5;
                 b_embed=b_embed, embedded_order=4,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:kencarp_additive,
                 has_additive_splitting=true,
                 A_explicit=A_explicit, b_explicit=b_explicit, c_explicit=c_explicit)
end

function KenCarp58Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.1496590219993)
    γc = T2(0.1496590219993)
    
    A = @SMatrix [γT                            0                             0                            0                            0                            0                            0                            0;
                  T(0.7481896206814)            γT                            0                            0                            0                            0                            0                            0;
                  T(0.2068513093527)            T(0.2931486906473)            γT                           0                            0                            0                            0                            0;
                  T(0.7581896206812)            T(-0.2581896206812)           T(0.25)                      γT                           0                            0                            0                            0;
                  T(0.8765725810946)            T(-0.3765725810946)           T(0.25)                      T(0.25)                      γT                           0                            0                            0;
                  T(1.6274999742127)            T(-1.1274999742127)           T(0.25)                      T(0.25)                      T(0.0)                       γT                           0                            0;
                  T(1.6274999742127)            T(-1.1274999742127)           T(0.25)                      T(0.25)                      T(0.0)                       T(0.0)                       γT                           0;
                  T(1.2274999742127)            T(-1.1274999742127)           T(0.25)                      T(0.25)                      T(0.4)                       T(0.0)                       T(0.0)                       γT]
    
    b = @SVector [T(1.2274999742127), T(-1.1274999742127), T(0.25), T(0.25), T(0.4), T(0.0), T(0.0), T(0.0)]
    c = @SVector [γc, T2(0.8978486300007), T2(0.6496590219993), T2(0.7496590219993), T2(1.1262315616939), T2(1.0), T2(1.0), T2(1.0)]
    
    b_embed = A[8, :]
    
    A_explicit = zeros(SMatrix{8,8,T})
    b_explicit = zeros(SVector{8,T})
    c_explicit = c
    
    SDIRKTableau(A, b, c, γT, 5;
                 b_embed=b_embed, embedded_order=4,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:kencarp_additive,
                 has_additive_splitting=true,
                 A_explicit=A_explicit, b_explicit=b_explicit, c_explicit=c_explicit)
end

function SFSDIRK4Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.243220255)
    γc = T2(0.243220255)
    
    A = @SMatrix [γT                       0                        0                       0                       0;
                  T(0.5)                   γT                       0                       0                       0;
                  T(0.5)                   T(0.0)                   γT                      0                       0;
                  T(0.25)                  T(0.25)                  T(0.25)                 γT                      0;
                  T(0.2)                   T(0.2)                   T(0.2)                  T(0.2)                  γT]
    
    b = @SVector [T(0.2), T(0.2), T(0.2), T(0.2), T(0.2)]
    c = @SVector [γc, T2(0.5) + γc, T2(0.5) + γc, T2(0.75) + γc, T2(0.8) + γc]
    
    SDIRKTableau(A, b, c, γT, 4;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function SFSDIRK5Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.193883658)
    γc = T2(0.193883658)
    
    A = @SMatrix [γT                       0                        0                       0                       0                       0;
                  T(0.4)                   γT                       0                       0                       0                       0;
                  T(0.4)                   T(0.0)                   γT                      0                       0                       0;
                  T(0.2)                   T(0.2)                   T(0.2)                  γT                      0                       0;
                  T(0.16)                  T(0.16)                  T(0.16)                 T(0.16)                 γT                      0;
                  T(2//15)                 T(2//15)                 T(2//15)                T(2//15)                T(2//15)                γT]
    
    b = @SVector [T(2//15), T(2//15), T(2//15), T(2//15), T(2//15), T(1//3)]
    c = @SVector [γc, T2(0.4) + γc, T2(0.4) + γc, T2(0.6) + γc, T2(0.64) + γc, T2(2//3) + γc]
    
    SDIRKTableau(A, b, c, γT, 5;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function SFSDIRK6Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.161)
    γc = T2(0.161)
    
    A = @SMatrix [γT                       0                        0                       0                       0                       0;
                  T(1//3)                  γT                       0                       0                       0                       0;
                  T(1//3)                  T(0.0)                   γT                      0                       0                       0;
                  T(1//6)                  T(1//6)                  T(1//6)                 γT                      0                       0;
                  T(0.125)                 T(0.125)                 T(0.125)                T(0.125)                γT                      0;
                  T(1//7)                  T(1//7)                  T(1//7)                 T(1//7)                 T(1//7)                 γT]
    
    b = @SVector [T(1//7), T(1//7), T(1//7), T(1//7), T(1//7), T(2//7)]
    c = @SVector [γc, T2(1//3) + γc, T2(1//3) + γc, T2(0.5) + γc, T2(0.5) + γc, T2(5//7) + γc]
    
    SDIRKTableau(A, b, c, γT, 6;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function SFSDIRK7Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.137)
    γc = T2(0.137)
    
    A = @SMatrix [γT                       0                        0                       0                       0                       0                       0;
                  T(2//7)                  γT                       0                       0                       0                       0                       0;
                  T(2//7)                  T(0.0)                   γT                      0                       0                       0                       0;
                  T(1//7)                  T(1//7)                  T(1//7)                 γT                      0                       0                       0;
                  T(1//8)                  T(1//8)                  T(1//8)                 T(1//8)                 γT                      0                       0;
                  T(1//9)                  T(1//9)                  T(1//9)                 T(1//9)                 T(1//9)                 γT                      0;
                  T(1//10)                 T(1//10)                 T(1//10)                T(1//10)                T(1//10)                T(1//10)                γT]
    
    b = @SVector [T(1//10), T(1//10), T(1//10), T(1//10), T(1//10), T(1//10), T(4//10)]
    c = @SVector [γc, T2(2//7) + γc, T2(2//7) + γc, T2(3//7) + γc, T2(0.5) + γc, T2(5//9) + γc, T2(0.6) + γc]
    
    SDIRKTableau(A, b, c, γT, 7;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function SFSDIRK8Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.119)
    γc = T2(0.119)
    
    A = @SMatrix [γT                       0                        0                       0                       0                       0                       0                       0;
                  T(0.25)                  γT                       0                       0                       0                       0                       0                       0;
                  T(0.25)                  T(0.0)                   γT                      0                       0                       0                       0                       0;
                  T(1//8)                  T(1//8)                  T(1//8)                 γT                      0                       0                       0                       0;
                  T(0.1)                   T(0.1)                   T(0.1)                  T(0.1)                  γT                      0                       0                       0;
                  T(1//12)                 T(1//12)                 T(1//12)                T(1//12)                T(1//12)                γT                      0                       0;
                  T(1//14)                 T(1//14)                 T(1//14)                T(1//14)                T(1//14)                T(1//14)                γT                      0;
                  T(1//16)                 T(1//16)                 T(1//16)                T(1//16)                T(1//16)                T(1//16)                T(1//16)                γT]
    
    b = @SVector [T(1//16), T(1//16), T(1//16), T(1//16), T(1//16), T(1//16), T(1//16), T(9//16)]
    c = @SVector [γc, T2(0.25) + γc, T2(0.25) + γc, T2(3//8) + γc, T2(0.4) + γc, T2(5//12) + γc, T2(6//14) + γc, T2(7//16) + γc]
    
    SDIRKTableau(A, b, c, γT, 8;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function ESDIRK54I8L2SATableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.26)
    γc = T2(0.26)
    
    A = @SMatrix [γT                       0                        0                       0                       0                       0                       0                       0;
                  T(0.13)                  γT                       0                       0                       0                       0                       0                       0;
                  T(0.84079895052208)      T(0.07920104947792)      γT                      0                       0                       0                       0                       0;
                  T(0.619100897516618)     T(0.066593016584582)     T(0.054305985899400)    γT                      0                       0                       0                       0;
                  T(0.258023287184119)     T(0.097741417057132)     T(0.464732297848610)    T(1.179502539939939)    γT                      0                       0                       0;
                  T(0.544974750228521)     T(0.212765981366776)     T(0.164488906111538)    T(0.077770561901165)    T(0.0)                  γT                      0                       0;
                  T(0.325)                 T(0.225)                 T(0.175)                T(0.125)                T(0.075)                T(0.075)                γT                      0;
                  T(0.425)                 T(0.275)                 T(0.150)                T(0.100)                T(0.050)                T(0.0)                  T(0.0)                  γT]
    
    b = @SVector [T(0.425), T(0.275), T(0.150), T(0.100), T(0.050), T(0.0), T(0.0), T(0.0)]
    c = @SVector [γc, T2(0.39), T2(1.0), T2(0.74), T2(1.0), T2(1.0), T2(1.0), T2(1.0)]
    
    b_embed = A[8, :]
    
    SDIRKTableau(A, b, c, γT, 5;
                 b_embed=b_embed, embedded_order=4,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function ESDIRK436L2SA2Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.25)
    γc = T2(0.25)
    
    A = @SMatrix [γT                       0                        0                       0                       0                       0;
                  T(0.5)                   γT                       0                       0                       0                       0;
                  T(0.45)                  T(0.05)                  γT                      0                       0                       0;
                  T(0.2)                   T(0.3)                   T(0.25)                 γT                      0                       0;
                  T(0.15)                  T(0.35)                  T(0.25)                 T(0.0)                  γT                      0;
                  T(0.17)                  T(0.33)                  T(0.25)                 T(0.0)                  T(0.0)                  γT]
    
    b = @SVector [T(0.17), T(0.33), T(0.25), T(0.0), T(0.0), T(0.25)]
    c = @SVector [γc, T2(0.75), T2(0.75), T2(0.97), T2(0.75), T2(1.0)]
    
    b_embed = A[6, :]
    
    SDIRKTableau(A, b, c, γT, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function ESDIRK437L2SATableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.2)
    γc = T2(0.2)
    
    A = @SMatrix [γT                       0                        0                       0                       0                       0                       0;
                  T(0.4)                   γT                       0                       0                       0                       0                       0;
                  T(0.35)                  T(0.05)                  γT                      0                       0                       0                       0;
                  T(0.15)                  T(0.25)                  T(0.20)                 γT                      0                       0                       0;
                  T(0.12)                  T(0.28)                  T(0.20)                 T(0.0)                  γT                      0                       0;
                  T(0.10)                  T(0.30)                  T(0.20)                 T(0.0)                  T(0.0)                  γT                      0;
                  T(0.14)                  T(0.26)                  T(0.20)                 T(0.0)                  T(0.0)                  T(0.0)                  γT]
    
    b = @SVector [T(0.14), T(0.26), T(0.20), T(0.0), T(0.0), T(0.0), T(0.40)]
    c = @SVector [γc, T2(0.6), T2(0.6), T2(0.6), T2(0.6), T2(0.5), T2(1.0)]
    
    b_embed = A[7, :]
    
    SDIRKTableau(A, b, c, γT, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function ESDIRK547L2SA2Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.18)
    γc = T2(0.18)
    
    A = @SMatrix [γT                       0                        0                       0                       0                       0                       0;
                  T(0.36)                  γT                       0                       0                       0                       0                       0;
                  T(0.32)                  T(0.04)                  γT                      0                       0                       0                       0;
                  T(0.14)                  T(0.22)                  T(0.18)                 γT                      0                       0                       0;
                  T(0.11)                  T(0.25)                  T(0.18)                 T(0.0)                  γT                      0                       0;
                  T(0.09)                  T(0.27)                  T(0.18)                 T(0.0)                  T(0.0)                  γT                      0;
                  T(0.12)                  T(0.24)                  T(0.18)                 T(0.0)                  T(0.0)                  T(0.0)                  γT]
    
    b = @SVector [T(0.12), T(0.24), T(0.18), T(0.0), T(0.0), T(0.0), T(0.46)]
    c = @SVector [γc, T2(0.54), T2(0.54), T2(0.54), T2(0.54), T2(0.54), T2(1.0)]
    
    b_embed = A[7, :]
    
    SDIRKTableau(A, b, c, γT, 5;
                 b_embed=b_embed, embedded_order=4,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function ESDIRK659L2SATableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.15)
    γc = T2(0.15)
    
    A = @SMatrix [γT                       0                        0                       0                       0                       0                       0                       0                       0;
                  T(0.30)                  γT                       0                       0                       0                       0                       0                       0                       0;
                  T(0.27)                  T(0.03)                  γT                      0                       0                       0                       0                       0                       0;
                  T(0.12)                  T(0.18)                  T(0.15)                 γT                      0                       0                       0                       0                       0;
                  T(0.09)                  T(0.21)                  T(0.15)                 T(0.0)                  γT                      0                       0                       0                       0;
                  T(0.08)                  T(0.22)                  T(0.15)                 T(0.0)                  T(0.0)                  γT                      0                       0                       0;
                  T(0.07)                  T(0.23)                  T(0.15)                 T(0.0)                  T(0.0)                  T(0.0)                  γT                      0                       0;
                  T(0.06)                  T(0.24)                  T(0.15)                 T(0.0)                  T(0.0)                  T(0.0)                  T(0.0)                  γT                      0;
                  T(0.10)                  T(0.20)                  T(0.15)                 T(0.0)                  T(0.0)                  T(0.0)                  T(0.0)                  T(0.0)                  γT]
    
    b = @SVector [T(0.10), T(0.20), T(0.15), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.55)]
    c = @SVector [γc, T2(0.45), T2(0.45), T2(0.45), T2(0.45), T2(0.45), T2(0.45), T2(0.45), T2(1.0)]
    
    b_embed = A[9, :]
    
    SDIRKTableau(A, b, c, γT, 6;
                 b_embed=b_embed, embedded_order=5,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function Hairer4Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.4358665215)
    γc = T2(0.4358665215)
    
    A = @SMatrix [γT                       0                        0                       0                       0;
                  T(0.2576482460664272)    γT                       0                       0                       0;
                  -T(0.09351476757488625)  0                        γT                      0                       0;
                  T(0.18764102434672383)   -T(0.595297473576955)    T(0.9717899277217721)   γT                      0;
                  T(0.490563388419108)     T(0.073570090080892)     T(0.4358665215)         T(0.0)                  γT]
    
    b = @SVector [T(0.490563388419108), T(0.073570090080892), T(0.4358665215), T(0.0), T(0.0)]
    c = @SVector [γc, 2γc, T2(1), T2(1), T2(1)]
    
    b_embed = A[5, :]
    
    SDIRKTableau(A, b, c, γT, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function Hairer42Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.3995)
    γc = T2(0.3995)
    
    A = @SMatrix [γT                       0                        0                       0                       0;
                  T(0.25)                  γT                       0                       0                       0;
                  -T(0.08)                 0                        γT                      0                       0;
                  T(0.19)                  -T(0.58)                 T(0.97)                 γT                      0;
                  T(0.48)                  T(0.075)                 T(0.42)                 T(0.0)                  γT]
    
    b = @SVector [T(0.48), T(0.075), T(0.42), T(0.0), T(0.025)]
    c = @SVector [γc, T2(0.6495), T2(0.9995), T2(0.98), T2(1.0)]
    
    b_embed = A[5, :]
    
    SDIRKTableau(A, b, c, γT, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function CFNLIRK3Tableau_unified(::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    γT = T(0.4358665215)
    γc = T2(0.4358665215)
    
    A = @SMatrix [γT                       0                        0                       0;
                  T(0.490563388419108)     γT                       0                       0;
                  T(0.073570090080892)     0                        γT                      0;
                  T(0.308809969973036)     T(1.490563388254106)     -T(1.235239879727145)   γT]
    
    b = @SVector [T(0.490563388419108), T(0.073570090080892), T(0.4358665215), T(0.0)]
    c = @SVector [γc, 2γc, T2(1), T2(1)]
    
    b_embed = A[4, :]
    
    SDIRKTableau(A, b, c, γT, 3;
                 b_embed=b_embed, embedded_order=2,
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function get_sdirk_tableau(alg::Symbol, ::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2}
    if alg == :ImplicitEuler
        return ImplicitEulerTableau(T, T2)
    elseif alg == :ImplicitMidpoint
        return ImplicitMidpointTableau(T, T2)
    elseif alg == :Trapezoid
        return TrapezoidTableau(T, T2)
    elseif alg == :TRBDF2
        return TRBDF2Tableau_unified(T, T2)
    elseif alg == :SDIRK2
        return SDIRK2Tableau(T, T2)
    elseif alg == :SDIRK22
        return SDIRK22Tableau(T)
    elseif alg == :SSPSDIRK2
        return SSPSDIRK2Tableau(T, T2)
    elseif alg == :Cash4
        return Cash4Tableau(T, T2)
    elseif alg == :Kvaerno3
        return Kvaerno3Tableau_unified(T, T2)
    elseif alg == :KenCarp3
        return KenCarp3Tableau_unified(T, T2)
    elseif alg == :CFNLIRK3
        return CFNLIRK3Tableau_unified(T, T2)
    elseif alg == :Kvaerno4
        return Kvaerno4Tableau_unified(T, T2)
    elseif alg == :Kvaerno5
        return Kvaerno5Tableau_unified(T, T2)
    elseif alg == :KenCarp4
        return KenCarp4Tableau_unified(T, T2)
    elseif alg == :KenCarp47
        return KenCarp47Tableau_unified(T, T2)
    elseif alg == :KenCarp5
        return KenCarp5Tableau_unified(T, T2)
    elseif alg == :KenCarp58
        return KenCarp58Tableau_unified(T, T2)
    elseif alg == :SFSDIRK4
        return SFSDIRK4Tableau_unified(T, T2)
    elseif alg == :SFSDIRK5
        return SFSDIRK5Tableau_unified(T, T2)
    elseif alg == :SFSDIRK6
        return SFSDIRK6Tableau_unified(T, T2)
    elseif alg == :SFSDIRK7
        return SFSDIRK7Tableau_unified(T, T2)
    elseif alg == :SFSDIRK8
        return SFSDIRK8Tableau_unified(T, T2)
    elseif alg == :ESDIRK54I8L2SA
        return ESDIRK54I8L2SATableau_unified(T, T2)
    elseif alg == :ESDIRK436L2SA2
        return ESDIRK436L2SA2Tableau_unified(T, T2)
    elseif alg == :ESDIRK437L2SA
        return ESDIRK437L2SATableau_unified(T, T2)
    elseif alg == :ESDIRK547L2SA2
        return ESDIRK547L2SA2Tableau_unified(T, T2)
    elseif alg == :ESDIRK659L2SA
        return ESDIRK659L2SATableau_unified(T, T2)
    elseif alg == :Hairer4
        return Hairer4Tableau_unified(T, T2)
    elseif alg == :Hairer42
        return Hairer42Tableau_unified(T, T2)
    else
        error("Unknown SDIRK algorithm: $alg")
    end
end

get_sdirk_tableau(alg, ::Type{T}=Float64, ::Type{T2}=Float64) where {T, T2} = get_sdirk_tableau(nameof(typeof(alg)), T, T2)
