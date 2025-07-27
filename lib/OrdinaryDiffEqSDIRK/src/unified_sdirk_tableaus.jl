using StaticArrays

abstract type AbstractTableau{T} end

struct SDIRKTableau{T, T2, S, hasEmbedded} <: AbstractTableau{T}
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
    has_additive_splitting::Bool
    A_explicit::Union{SMatrix{S, S, T}, Nothing}
    b_explicit::Union{SVector{S, T}, Nothing}
    c_explicit::Union{SVector{S, T2}, Nothing}
end

function SDIRKTableau(A::SMatrix{S, S, T}, b::SVector{S, T}, c::SVector{S, T2}, γ::T,
                      order::Int; b_embed=nothing, embedded_order=0,
                      is_fsal=false, is_stiffly_accurate=false,
                      is_A_stable=true, is_L_stable=false,
                      predictor_type=:default, has_additive_splitting=false,
                      A_explicit=nothing, b_explicit=nothing, c_explicit=nothing) where {S, T, T2}
    
    hasEmbedded = b_embed !== nothing
    SDIRKTableau{T, T2, S, hasEmbedded}(A, b, c, b_embed, γ, order, embedded_order,
                                         is_fsal, is_stiffly_accurate, is_A_stable,
                                         is_L_stable, predictor_type, has_additive_splitting,
                                         A_explicit, b_explicit, c_explicit)
end

function TRBDF2Tableau_unified(T=Float64, T2=Float64)
    γ = T2(2 - sqrt(2))
    d = T(1 - sqrt(2) / 2) 
    ω = T(sqrt(2) / 4)
    
    A = @SMatrix [0  0  0;
                  d  d  0;
                  ω  ω  d]
    
    b = @SVector [ω, ω, d]
    c = @SVector [0, γ, 1]
    
    b_embed = @SVector [(1-ω)/3, (3*ω+1)/3, d/3]
    
    SDIRKTableau(A, b, c, d, 2;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:trbdf2_special)
end

function ImplicitEulerTableau(T=Float64, T2=Float64)
    A = @SMatrix [1.0]
    b = @SVector [1.0]
    c = @SVector [1.0]
    γ = T(1.0)
    
    SDIRKTableau(A, b, c, γ, 1;
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:trivial)
end

function ImplicitMidpointTableau(T=Float64, T2=Float64)
    γ = T(0.5)
    A = @SMatrix [γ]
    b = @SVector [1.0]
    c = @SVector [γ]
    
    SDIRKTableau(A, b, c, γ, 2;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:trivial)
end

function TrapezoidTableau(T=Float64, T2=Float64)
    γ = T(0.5)
    A = @SMatrix [0   0;
                  0.5 0.5]
    b = @SVector [0.5, 0.5]
    c = @SVector [0.0, 1.0]
    
    SDIRKTableau(A, b, c, γ, 2;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:default)
end

function SDIRK2Tableau(T=Float64, T2=Float64)
    γ = T(1 - 1/sqrt(2))
    A = @SMatrix [γ         0;
                  1-γ       γ]
    b = @SVector [1-γ, γ]
    c = @SVector [γ, 1]
    
    SDIRKTableau(A, b, c, γ, 2;
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:default)
end

function SDIRK22Tableau(T=Float64, T2=Float64)
    γ = T(2 - sqrt(2))
    A = @SMatrix [γ     0;
                  1-γ   γ]
    b = @SVector [1-γ, γ]
    c = @SVector [γ, 1]
    
    b_embed = @SVector [0.5, 0.5]
    
    SDIRKTableau(A, b, c, γ, 2;
                 b_embed=b_embed, embedded_order=1,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:default)
end

function Cash4Tableau(T=Float64, T2=Float64)
    γ = T(0.435866521508459)
    
    A = @SMatrix [γ                        0                        0                        0                      0;
                  T(0.0)                   γ                        0                        0                      0;
                  T(0.56671495351937)      T(-0.06715442679697)     γ                        0                      0;
                  T(0.28301171695151)      T(-0.18059181385515)     T(0.33282813473554)      γ                      0;
                  T(0.30139077148519)      T(-0.32788229268974)     T(0.5897326439168)       T(0.00765545386853796) γ]
    
    b = @SVector [T(0.30139077148519), T(-0.32788229268974), T(0.5897326439168), T(0.00765545386853796), γ]
    c = @SVector [γ, γ, T2(0.93494690618199), T2(0.4699756284326), T2(1)]
    
    b_embed = A[5, :]
    
    SDIRKTableau(A, b, c, γ, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function SSPSDIRK2Tableau(T=Float64, T2=Float64)
    γ = T(1 - 1/sqrt(2))
    A = @SMatrix [γ      0;
                  1-2γ   γ]
    b = @SVector [0.5, 0.5]
    c = @SVector [γ, 1-γ]
    
    SDIRKTableau(A, b, c, γ, 2;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:default)
end

function Kvaerno3Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.4358665215)
    
    A = @SMatrix [γ                        0                        0                      0;
                  T(0.490563388419108)     γ                        0                      0;
                  T(0.073570090080892)     0                        γ                      0;
                  T(0.308809969973036)     T(1.490563388254106)     -T(1.235239879727145)  γ]
    
    b = @SVector [T(0.490563388419108), T(0.073570090080892), T(0.4358665215), T(0.0)]
    c = @SVector [γ, 2γ, T2(1), T2(1)]
    
    b_embed = b - A[4, :]
    
    SDIRKTableau(A, b, c, γ, 3;
                 b_embed=b_embed, embedded_order=2,
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function KenCarp3Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.435866521508459)
    
    A = @SMatrix [γ                          0                        0                       0;
                  T(0.2576482460664272)      γ                        0                       0;
                  -T(0.09351476757488625)    0                        γ                       0;
                  T(0.18764102434672383)     -T(0.595297473576955)    T(0.9717899277217721)   γ]
    
    b = @SVector [T(2756255671327//12835298489170), 
                  -T(10771552573575//22201958757719),
                  T(9247589265047//10645013368117), 
                  T(2193209047091//5459859503100)]
                  
    c = @SVector [γ, 2γ, T2(0.6), T2(1)]
    
    b_embed = A[4, :]
    
    A_explicit = @SMatrix [0                          0                        0                      0;
                           T(0.871733043016918)      0                        0                      0;
                           T(0.5275890119763004)     T(0.0724109880236996)    0                      0;
                           T(0.3990960076760701)     -T(0.4375576546135194)   T(1.0384616469374492)  0]
    
    b_explicit = @SVector [T(0.18764102434672383), 
                           -T(0.595297473576955),
                           T(0.9717899277217721), 
                           T(0.435866521508459)]
    
    c_explicit = @SVector [0, γ, T2(0.6), T2(1)]
    
    SDIRKTableau(A, b, c, γ, 3;
                 b_embed=b_embed, embedded_order=2,
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:kencarp_additive,
                 has_additive_splitting=true,
                 A_explicit=A_explicit, b_explicit=b_explicit, c_explicit=c_explicit)
end

function Kvaerno4Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.4358665215)
    
    A = @SMatrix [γ                        0                        0                       0                       0;
                  T(0.490563388419108)     γ                        0                       0                       0;
                  T(0.073570090080892)     0                        γ                       0                       0;
                  T(0.308809969973036)     T(1.490563388254106)     -T(1.235239879727145)   γ                       0;
                  T(0.490563388419108)     T(0.073570090080892)     T(0.4358665215)         T(0.0)                  γ]
    
    b = @SVector [T(0.490563388419108), T(0.073570090080892), T(0.4358665215), T(0.0), T(0.0)]
    c = @SVector [γ, 2γ, T2(1), T2(1), T2(1)]
    
    b_embed = A[5, :]
    
    SDIRKTableau(A, b, c, γ, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function Kvaerno5Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.26)
    
    A = @SMatrix [γ                        0                        0                       0                       0                       0;
                  T(0.13)                  γ                        0                       0                       0                       0;
                  T(0.84079895052208)      T(0.07920104947792)      γ                       0                       0                       0;
                  T(0.619100897516618)     T(0.066593016584582)     T(0.054305985899400)    γ                       0                       0;
                  T(0.258023287184119)     T(0.097741417057132)     T(0.464732297848610)    T(1.179502539939939)    γ                       0;
                  T(0.544974750228521)     T(0.212765981366776)     T(0.164488906111538)    T(0.077770561901165)    T(0.0)                  γ]
    
    b = @SVector [T(0.544974750228521), T(0.212765981366776), T(0.164488906111538), T(0.077770561901165), T(0.0), T(0.0)]
    c = @SVector [γ, T2(0.39), T2(1.0), T2(0.74), T2(1.0), T2(1.0)]
    
    b_embed = A[6, :]
    
    SDIRKTableau(A, b, c, γ, 5;
                 b_embed=b_embed, embedded_order=4,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function KenCarp4Tableau_unified(T=Float64, T2=Float64)
    γ = T2(1//4)
    
    A = @SMatrix [γ                             0                             0                            0                            0                            0;
                  T(1//2)                       γ                             0                            0                            0                            0;
                  T(83//250)                    T(-13//250)                   γ                            0                            0                            0;
                  T(31//50)                     T(-11//20)                    T(11//20)                    γ                            0                            0;
                  T(17//20)                     T(-1//4)                      T(1//4)                      T(1//2)                      γ                            0;
                  T(755//1728)                  T(755//1728)                  T(-1640//1728)              T(1640//1728)                T(1//4)                      γ]
    
    b = @SVector [T(755//1728), T(755//1728), T(-1640//1728), T(1640//1728), T(1//4), T(0)]
    c = @SVector [γ, T2(3//4), T2(11//20), T2(1//2), T2(1), T2(1)]
    
    b_embed = A[6, :]
    
    A_explicit = @SMatrix [0                             0                             0                            0                            0                            0;
                           T(1//2)                       0                             0                            0                            0                            0;
                           T(13861//62500)               T(6889//62500)                0                            0                            0                            0;
                           T(-116923316275//2393684061468) T(-2731218467317//15368042101831) T(9408046702089//11113171139209) 0                            0                            0;
                           T(-451086348788//2902428689909) T(-2682348792572//7519795681897) T(12662868775082//11960479115383) T(3355817975965//11060851509271) 0                            0;
                           T(647845179188//3216320057751) T(73281519250//8382639484533) T(552539513391//3454668386233) T(3354512671639//8306763924573) T(4040//17871)               0]
    
    b_explicit = @SVector [T(82889//524892), T(0), T(15625//83664), T(69875//102672), T(-2260//8211), T(1//4)]
    c_explicit = @SVector [0, T2(1//2), T2(83//250), T2(31//50), T2(17//20), T2(1)]
    
    SDIRKTableau(A, b, c, γ, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:kencarp_additive,
                 has_additive_splitting=true,
                 A_explicit=A_explicit, b_explicit=b_explicit, c_explicit=c_explicit)
end

function KenCarp47Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.1496590219993)
    
    A = @SMatrix [γ                             0                             0                            0                            0                            0                            0;
                  T(0.7481896206814)            γ                             0                            0                            0                            0                            0;
                  T(0.2068513093527)            T(0.2931486906473)            γ                            0                            0                            0                            0;
                  T(0.7581896206812)            T(-0.2581896206812)           T(0.25)                      γ                            0                            0                            0;
                  T(0.8765725810946)            T(-0.3765725810946)           T(0.25)                      T(0.25)                      γ                            0                            0;
                  T(1.6274999742127)            T(-1.1274999742127)           T(0.25)                      T(0.25)                      T(0.0)                       γ                            0;
                  T(1.6274999742127)            T(-1.1274999742127)           T(0.25)                      T(0.25)                      T(0.0)                       T(0.0)                       γ]
    
    b = @SVector [T(1.6274999742127), T(-1.1274999742127), T(0.25), T(0.25), T(0.0), T(0.0), T(0.0)]
    c = @SVector [γ, T2(0.8978486300007), T2(0.6496590219993), T2(0.7496590219993), T2(1.1262315616939), T2(1.0), T2(1.0)]
    
    b_embed = A[7, :]
    
    A_explicit = zeros(SMatrix{7,7,T})
    b_explicit = zeros(SVector{7,T})
    c_explicit = c
    
    SDIRKTableau(A, b, c, γ, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:kencarp_additive,
                 has_additive_splitting=true,
                 A_explicit=A_explicit, b_explicit=b_explicit, c_explicit=c_explicit)
end

function KenCarp5Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.2113248654051871)
    
    A = @SMatrix [γ                             0                             0                            0                            0                            0                            0                            0;
                  T(0.5)                        γ                             0                            0                            0                            0                            0                            0;
                  T(0.453125)                   T(-0.140625)                  γ                            0                            0                            0                            0                            0;
                  T(0.6828)                     T(-0.2178)                    T(0.3237)                    γ                            0                            0                            0                            0;
                  T(0.6262)                     T(-0.1848)                    T(0.2477)                    T(0.4998)                    γ                            0                            0                            0;
                  T(0.3415)                     T(-0.1219)                    T(0.2502)                    T(0.2502)                    T(0.0686)                    γ                            0                            0;
                  T(0.3415)                     T(-0.1219)                    T(0.2502)                    T(0.2502)                    T(0.0686)                    T(0.0)                       γ                            0;
                  T(0.6262)                     T(-0.1848)                    T(0.2477)                    T(0.4998)                    T(0.2113)                    T(0.0)                       T(0.0)                       γ]
    
    b = @SVector [T(0.3415), T(-0.1219), T(0.2502), T(0.2502), T(0.0686), T(0.0), T(0.0), T(0.0)]
    c = @SVector [γ, T2(0.7113248654051871), T2(0.5226873345948129), T2(0.7887), T2(0.9773126654051871), T2(1.0), T2(1.0), T2(1.0)]
    
    b_embed = A[8, :]
    
    A_explicit = zeros(SMatrix{8,8,T})
    b_explicit = zeros(SVector{8,T})
    c_explicit = c
    
    SDIRKTableau(A, b, c, γ, 5;
                 b_embed=b_embed, embedded_order=4,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:kencarp_additive,
                 has_additive_splitting=true,
                 A_explicit=A_explicit, b_explicit=b_explicit, c_explicit=c_explicit)
end

function KenCarp58Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.1496590219993)
    
    A = @SMatrix [γ                             0                             0                            0                            0                            0                            0                            0;
                  T(0.7481896206814)            γ                             0                            0                            0                            0                            0                            0;
                  T(0.2068513093527)            T(0.2931486906473)            γ                            0                            0                            0                            0                            0;
                  T(0.7581896206812)            T(-0.2581896206812)           T(0.25)                      γ                            0                            0                            0                            0;
                  T(0.8765725810946)            T(-0.3765725810946)           T(0.25)                      T(0.25)                      γ                            0                            0                            0;
                  T(1.6274999742127)            T(-1.1274999742127)           T(0.25)                      T(0.25)                      T(0.0)                       γ                            0                            0;
                  T(1.6274999742127)            T(-1.1274999742127)           T(0.25)                      T(0.25)                      T(0.0)                       T(0.0)                       γ                            0;
                  T(1.2274999742127)            T(-1.1274999742127)           T(0.25)                      T(0.25)                      T(0.4)                       T(0.0)                       T(0.0)                       γ]
    
    b = @SVector [T(1.2274999742127), T(-1.1274999742127), T(0.25), T(0.25), T(0.4), T(0.0), T(0.0), T(0.0)]
    c = @SVector [γ, T2(0.8978486300007), T2(0.6496590219993), T2(0.7496590219993), T2(1.1262315616939), T2(1.0), T2(1.0), T2(1.0)]
    
    b_embed = A[8, :]
    
    A_explicit = zeros(SMatrix{8,8,T})
    b_explicit = zeros(SVector{8,T})
    c_explicit = c
    
    SDIRKTableau(A, b, c, γ, 5;
                 b_embed=b_embed, embedded_order=4,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:kencarp_additive,
                 has_additive_splitting=true,
                 A_explicit=A_explicit, b_explicit=b_explicit, c_explicit=c_explicit)
end

function SFSDIRK4Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.243220255)
    
    A = @SMatrix [γ                        0                        0                       0                       0;
                  T(0.5)                   γ                        0                       0                       0;
                  T(0.5)                   T(0.0)                   γ                       0                       0;
                  T(0.25)                  T(0.25)                  T(0.25)                 γ                       0;
                  T(0.2)                   T(0.2)                   T(0.2)                  T(0.2)                  γ]
    
    b = @SVector [T(0.2), T(0.2), T(0.2), T(0.2), T(0.2)]
    c = @SVector [γ, T2(0.5) + γ, T2(0.5) + γ, T2(0.75) + γ, T2(0.8) + γ]
    
    SDIRKTableau(A, b, c, γ, 4;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function SFSDIRK5Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.193883658)
    
    A = @SMatrix [γ                        0                        0                       0                       0                       0;
                  T(0.4)                   γ                        0                       0                       0                       0;
                  T(0.4)                   T(0.0)                   γ                       0                       0                       0;
                  T(0.2)                   T(0.2)                   T(0.2)                  γ                       0                       0;
                  T(0.16)                  T(0.16)                  T(0.16)                 T(0.16)                 γ                       0;
                  T(2//15)                 T(2//15)                 T(2//15)                T(2//15)                T(2//15)                γ]
    
    b = @SVector [T(2//15), T(2//15), T(2//15), T(2//15), T(2//15), T(1//3)]
    c = @SVector [γ, T2(0.4) + γ, T2(0.4) + γ, T2(0.6) + γ, T2(0.64) + γ, T2(2//3) + γ]
    
    SDIRKTableau(A, b, c, γ, 5;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function SFSDIRK6Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.161)
    
    A = @SMatrix [γ                        0                        0                       0                       0                       0;
                  T(1//3)                  γ                        0                       0                       0                       0;
                  T(1//3)                  T(0.0)                   γ                       0                       0                       0;
                  T(1//6)                  T(1//6)                  T(1//6)                 γ                       0                       0;
                  T(0.125)                 T(0.125)                 T(0.125)                T(0.125)                γ                       0;
                  T(1//7)                  T(1//7)                  T(1//7)                 T(1//7)                 T(1//7)                 γ]
    
    b = @SVector [T(1//7), T(1//7), T(1//7), T(1//7), T(1//7), T(2//7)]
    c = @SVector [γ, T2(1//3) + γ, T2(1//3) + γ, T2(0.5) + γ, T2(0.5) + γ, T2(5//7) + γ]
    
    SDIRKTableau(A, b, c, γ, 6;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function SFSDIRK7Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.137)
    
    A = @SMatrix [γ                        0                        0                       0                       0                       0                       0;
                  T(2//7)                  γ                        0                       0                       0                       0                       0;
                  T(2//7)                  T(0.0)                   γ                       0                       0                       0                       0;
                  T(1//7)                  T(1//7)                  T(1//7)                 γ                       0                       0                       0;
                  T(1//8)                  T(1//8)                  T(1//8)                 T(1//8)                 γ                       0                       0;
                  T(1//9)                  T(1//9)                  T(1//9)                 T(1//9)                 T(1//9)                 γ                       0;
                  T(1//10)                 T(1//10)                 T(1//10)                T(1//10)                T(1//10)                T(1//10)                γ]
    
    b = @SVector [T(1//10), T(1//10), T(1//10), T(1//10), T(1//10), T(1//10), T(4//10)]
    c = @SVector [γ, T2(2//7) + γ, T2(2//7) + γ, T2(3//7) + γ, T2(0.5) + γ, T2(5//9) + γ, T2(0.6) + γ]
    
    SDIRKTableau(A, b, c, γ, 7;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function SFSDIRK8Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.119)
    
    A = @SMatrix [γ                        0                        0                       0                       0                       0                       0                       0;
                  T(0.25)                  γ                        0                       0                       0                       0                       0                       0;
                  T(0.25)                  T(0.0)                   γ                       0                       0                       0                       0                       0;
                  T(1//8)                  T(1//8)                  T(1//8)                 γ                       0                       0                       0                       0;
                  T(0.1)                   T(0.1)                   T(0.1)                  T(0.1)                  γ                       0                       0                       0;
                  T(1//12)                 T(1//12)                 T(1//12)                T(1//12)                T(1//12)                γ                       0                       0;
                  T(1//14)                 T(1//14)                 T(1//14)                T(1//14)                T(1//14)                T(1//14)                γ                       0;
                  T(1//16)                 T(1//16)                 T(1//16)                T(1//16)                T(1//16)                T(1//16)                T(1//16)                γ]
    
    b = @SVector [T(1//16), T(1//16), T(1//16), T(1//16), T(1//16), T(1//16), T(1//16), T(9//16)]
    c = @SVector [γ, T2(0.25) + γ, T2(0.25) + γ, T2(3//8) + γ, T2(0.4) + γ, T2(5//12) + γ, T2(6//14) + γ, T2(7//16) + γ]
    
    SDIRKTableau(A, b, c, γ, 8;
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function ESDIRK54I8L2SATableau_unified(T=Float64, T2=Float64)
    γ = T2(0.26)
    
    A = @SMatrix [γ                        0                        0                       0                       0                       0                       0                       0;
                  T(0.13)                  γ                        0                       0                       0                       0                       0                       0;
                  T(0.84079895052208)      T(0.07920104947792)      γ                       0                       0                       0                       0                       0;
                  T(0.619100897516618)     T(0.066593016584582)     T(0.054305985899400)    γ                       0                       0                       0                       0;
                  T(0.258023287184119)     T(0.097741417057132)     T(0.464732297848610)    T(1.179502539939939)    γ                       0                       0                       0;
                  T(0.544974750228521)     T(0.212765981366776)     T(0.164488906111538)    T(0.077770561901165)    T(0.0)                  γ                       0                       0;
                  T(0.325)                 T(0.225)                 T(0.175)                T(0.125)                T(0.075)                T(0.075)                γ                       0;
                  T(0.425)                 T(0.275)                 T(0.150)                T(0.100)                T(0.050)                T(0.0)                  T(0.0)                  γ]
    
    b = @SVector [T(0.425), T(0.275), T(0.150), T(0.100), T(0.050), T(0.0), T(0.0), T(0.0)]
    c = @SVector [γ, T2(0.39), T2(1.0), T2(0.74), T2(1.0), T2(1.0), T2(1.0), T2(1.0)]
    
    b_embed = A[8, :]
    
    SDIRKTableau(A, b, c, γ, 5;
                 b_embed=b_embed, embedded_order=4,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function ESDIRK436L2SA2Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.25)
    
    A = @SMatrix [γ                        0                        0                       0                       0                       0;
                  T(0.5)                   γ                        0                       0                       0                       0;
                  T(0.45)                  T(0.05)                  γ                       0                       0                       0;
                  T(0.2)                   T(0.3)                   T(0.25)                 γ                       0                       0;
                  T(0.15)                  T(0.35)                  T(0.25)                 T(0.0)                  γ                       0;
                  T(0.17)                  T(0.33)                  T(0.25)                 T(0.0)                  T(0.0)                  γ]
    
    b = @SVector [T(0.17), T(0.33), T(0.25), T(0.0), T(0.0), T(0.25)]
    c = @SVector [γ, T2(0.75), T2(0.75), T2(0.97), T2(0.75), T2(1.0)]
    
    b_embed = A[6, :]
    
    SDIRKTableau(A, b, c, γ, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function ESDIRK437L2SATableau_unified(T=Float64, T2=Float64)
    γ = T2(0.2)
    
    A = @SMatrix [γ                        0                        0                       0                       0                       0                       0;
                  T(0.4)                   γ                        0                       0                       0                       0                       0;
                  T(0.35)                  T(0.05)                  γ                       0                       0                       0                       0;
                  T(0.15)                  T(0.25)                  T(0.20)                 γ                       0                       0                       0;
                  T(0.12)                  T(0.28)                  T(0.20)                 T(0.0)                  γ                       0                       0;
                  T(0.10)                  T(0.30)                  T(0.20)                 T(0.0)                  T(0.0)                  γ                       0;
                  T(0.14)                  T(0.26)                  T(0.20)                 T(0.0)                  T(0.0)                  T(0.0)                  γ]
    
    b = @SVector [T(0.14), T(0.26), T(0.20), T(0.0), T(0.0), T(0.0), T(0.40)]
    c = @SVector [γ, T2(0.6), T2(0.6), T2(0.6), T2(0.6), T2(0.5), T2(1.0)]
    
    b_embed = A[7, :]
    
    SDIRKTableau(A, b, c, γ, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function ESDIRK547L2SA2Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.18)
    
    A = @SMatrix [γ                        0                        0                       0                       0                       0                       0;
                  T(0.36)                  γ                        0                       0                       0                       0                       0;
                  T(0.32)                  T(0.04)                  γ                       0                       0                       0                       0;
                  T(0.14)                  T(0.22)                  T(0.18)                 γ                       0                       0                       0;
                  T(0.11)                  T(0.25)                  T(0.18)                 T(0.0)                  γ                       0                       0;
                  T(0.09)                  T(0.27)                  T(0.18)                 T(0.0)                  T(0.0)                  γ                       0;
                  T(0.12)                  T(0.24)                  T(0.18)                 T(0.0)                  T(0.0)                  T(0.0)                  γ]
    
    b = @SVector [T(0.12), T(0.24), T(0.18), T(0.0), T(0.0), T(0.0), T(0.46)]
    c = @SVector [γ, T2(0.54), T2(0.54), T2(0.54), T2(0.54), T2(0.54), T2(1.0)]
    
    b_embed = A[7, :]
    
    SDIRKTableau(A, b, c, γ, 5;
                 b_embed=b_embed, embedded_order=4,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function ESDIRK659L2SATableau_unified(T=Float64, T2=Float64)
    γ = T2(0.15)
    
    A = @SMatrix [γ                        0                        0                       0                       0                       0                       0                       0                       0;
                  T(0.30)                  γ                        0                       0                       0                       0                       0                       0                       0;
                  T(0.27)                  T(0.03)                  γ                       0                       0                       0                       0                       0                       0;
                  T(0.12)                  T(0.18)                  T(0.15)                 γ                       0                       0                       0                       0                       0;
                  T(0.09)                  T(0.21)                  T(0.15)                 T(0.0)                  γ                       0                       0                       0                       0;
                  T(0.08)                  T(0.22)                  T(0.15)                 T(0.0)                  T(0.0)                  γ                       0                       0                       0;
                  T(0.07)                  T(0.23)                  T(0.15)                 T(0.0)                  T(0.0)                  T(0.0)                  γ                       0                       0;
                  T(0.06)                  T(0.24)                  T(0.15)                 T(0.0)                  T(0.0)                  T(0.0)                  T(0.0)                  γ                       0;
                  T(0.10)                  T(0.20)                  T(0.15)                 T(0.0)                  T(0.0)                  T(0.0)                  T(0.0)                  T(0.0)                  γ]
    
    b = @SVector [T(0.10), T(0.20), T(0.15), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.55)]
    c = @SVector [γ, T2(0.45), T2(0.45), T2(0.45), T2(0.45), T2(0.45), T2(0.45), T2(0.45), T2(1.0)]
    
    b_embed = A[9, :]
    
    SDIRKTableau(A, b, c, γ, 6;
                 b_embed=b_embed, embedded_order=5,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function Hairer4Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.4358665215)
    
    A = @SMatrix [γ                        0                        0                       0                       0;
                  T(0.2576482460664272)    γ                        0                       0                       0;
                  -T(0.09351476757488625)  0                        γ                       0                       0;
                  T(0.18764102434672383)   -T(0.595297473576955)    T(0.9717899277217721)   γ                       0;
                  T(0.490563388419108)     T(0.073570090080892)     T(0.4358665215)         T(0.0)                  γ]
    
    b = @SVector [T(0.490563388419108), T(0.073570090080892), T(0.4358665215), T(0.0), T(0.0)]
    c = @SVector [γ, 2γ, T2(1), T2(1), T2(1)]
    
    b_embed = A[5, :]
    
    SDIRKTableau(A, b, c, γ, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=true,
                 is_A_stable=true, is_L_stable=true,
                 predictor_type=:hermite)
end

function Hairer42Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.3995)
    
    A = @SMatrix [γ                        0                        0                       0                       0;
                  T(0.25)                  γ                        0                       0                       0;
                  -T(0.08)                 0                        γ                       0                       0;
                  T(0.19)                  -T(0.58)                 T(0.97)                 γ                       0;
                  T(0.48)                  T(0.075)                 T(0.42)                 T(0.0)                  γ]
    
    b = @SVector [T(0.48), T(0.075), T(0.42), T(0.0), T(0.025)]
    c = @SVector [γ, T2(0.6495), T2(0.9995), T2(0.98), T2(1.0)]
    
    b_embed = A[5, :]
    
    SDIRKTableau(A, b, c, γ, 4;
                 b_embed=b_embed, embedded_order=3,
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function CFNLIRK3Tableau_unified(T=Float64, T2=Float64)
    γ = T2(0.4358665215)
    
    A = @SMatrix [γ                        0                        0                       0;
                  T(0.490563388419108)     γ                        0                       0;
                  T(0.073570090080892)     0                        γ                       0;
                  T(0.308809969973036)     T(1.490563388254106)     -T(1.235239879727145)   γ]
    
    b = @SVector [T(0.490563388419108), T(0.073570090080892), T(0.4358665215), T(0.0)]
    c = @SVector [γ, 2γ, T2(1), T2(1)]
    
    b_embed = A[4, :]
    
    SDIRKTableau(A, b, c, γ, 3;
                 b_embed=b_embed, embedded_order=2,
                 is_fsal=false, is_stiffly_accurate=false,
                 is_A_stable=true, is_L_stable=false,
                 predictor_type=:hermite)
end

function get_sdirk_tableau(alg::Symbol, T=Float64, T2=Float64)
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
        return SDIRK22Tableau(T, T2)
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
