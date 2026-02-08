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

function SDIRKTableau(
        A::SMatrix{S, S, T}, b::SVector{S, T}, c::SVector{S, T2}, γ::T,
        order::Int; b_embed = nothing, embedded_order = 0,
        is_fsal = false, is_stiffly_accurate = false,
        is_A_stable = true, is_L_stable = false,
        predictor_type = :default, has_additive_splitting = false,
        A_explicit = nothing, b_explicit = nothing, c_explicit = nothing,
        α_pred = nothing, has_spice_error = false
    ) where {S, T, T2}

    hasEmbedded = b_embed !== nothing
    hasAdditiveSplitting = has_additive_splitting
    hasExplicit = A_explicit !== nothing
    hasPred = α_pred !== nothing
    return SDIRKTableau{T, T2, S, hasEmbedded, hasAdditiveSplitting, hasExplicit, hasPred}(
        A, b, c, b_embed, γ, order, embedded_order,
        is_fsal, is_stiffly_accurate, is_A_stable,
        is_L_stable, predictor_type,
        A_explicit, b_explicit, c_explicit, α_pred, has_spice_error
    )
end

function TRBDF2Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    γ = T2(2 - sqrt(2))
    d = T(1 - sqrt(2) / 2)
    ω = T(sqrt(2) / 4)

    A = @SMatrix [
        0  0  0;
        d  d  0;
        ω  ω  d
    ]

    b = @SVector [ω, ω, d]
    c = @SVector [0, γ, 1]

    # btilde = bhat - b, NOT bhat itself
    btilde1 = T((1 - sqrt(2)) / 3)
    btilde2 = T(1 // 3)
    btilde3 = T((sqrt(2) - 2) / 3)
    b_embed = @SVector [btilde1, btilde2, btilde3]

    α1 = T2(-sqrt(2) / 2)
    α2 = T2(1 + sqrt(2) / 2)
    α_pred = @SMatrix T2[
        0 0 0;
        0 0 0;
        α1 α2 0
    ]

    return SDIRKTableau(
        A, b, c, d, 2;
        b_embed = b_embed, embedded_order = 3,
        is_fsal = false, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :trbdf2_special,
        α_pred = α_pred
    )
end

function ImplicitEulerTableau(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    A = @SMatrix [T(1.0)]
    b = @SVector [T(1.0)]
    c = @SVector [T2(1.0)]
    γ = T(1.0)

    return SDIRKTableau(
        A, b, c, γ, 1;
        is_fsal = false, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :trivial,
        has_spice_error = true
    )
end

function ImplicitMidpointTableau(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    γT = T(0.5)
    γc = T2(0.5)
    A = @SMatrix [γT]
    b = @SVector [T(1.0)]
    c = @SVector [γc]

    return SDIRKTableau(
        A, b, c, γT, 2;
        is_fsal = false, is_stiffly_accurate = false,
        is_A_stable = true, is_L_stable = false,
        predictor_type = :trivial
    )
end

function TrapezoidTableau(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    γT = T(0.5)
    A = @SMatrix [
        T(0)   T(0);
        T(0.5) T(0.5)
    ]
    b = @SVector [T(0.5), T(0.5)]
    c = @SVector [T2(0.0), T2(1.0)]

    return SDIRKTableau(
        A, b, c, γT, 2;
        is_fsal = false, is_stiffly_accurate = false,
        is_A_stable = true, is_L_stable = false,
        predictor_type = :default,
        has_spice_error = true
    )
end

function SDIRK2Tableau(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    γT = T(1 - 1 / sqrt(2))
    γc = T2(1 - 1 / sqrt(2))
    # Master z-space: stage 1 tmp=uprev, stage 2 tmp=uprev-z₁
    # So A[2,1] = -1 in z-space
    A = @SMatrix [
        γT     T(0);
        T(-1)  γT
    ]
    # u = uprev + z₁/2 + z₂/2
    b = @SVector [T(1 // 2), T(1 // 2)]
    c = @SVector [γc, γc]  # master uses same c for both stages
    # btilde = z₁/2 - z₂/2
    b_embed = @SVector [T(1 // 2), T(-1 // 2)]

    return SDIRKTableau(
        A, b, c, γT, 2;
        b_embed = b_embed, embedded_order = 1,
        is_fsal = false, is_stiffly_accurate = false,
        is_A_stable = true, is_L_stable = false,
        predictor_type = :default
    )
end

function SSPSDIRK2Tableau(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    # Master hardcodes γ = 1/4
    γT = T(1 // 4)
    γc = T2(1 // 4)
    # Stage 1: tmp = uprev, c = 1
    # Stage 2: tmp = uprev + z₁/2, c = 1
    # u = tmp + z₂/2 = uprev + z₁/2 + z₂/2
    A = @SMatrix [
        γT        T(0);
        T(1 // 2)   γT
    ]
    b = @SVector [T(1 // 2), T(1 // 2)]
    c = @SVector [T2(1), T2(1)]

    return SDIRKTableau(
        A, b, c, γT, 2;
        is_fsal = false, is_stiffly_accurate = false,
        is_A_stable = true, is_L_stable = false,
        predictor_type = :default
    )
end

function Cash4Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    γT = T(0.435866521508)
    γc = T2(0.435866521508)
    a21 = T(-1.1358665215)
    a31 = T(1.08543330679)
    a32 = -T(0.721299828287)
    a41 = T(0.416349501547)
    a42 = T(0.190984004184)
    a43 = -T(0.118643265417)
    a51 = T(0.896869652944)
    a52 = T(0.0182725272734)
    a53 = -T(0.0845900310706)
    a54 = -T(0.266418670647)

    b1hat2 = T(0.77669193291)
    b2hat2 = T(0.0297472791484)
    b3hat2 = -T(0.0267440239074)
    b4hat2 = T(0.220304811849)

    A = @SMatrix [
        γT    T(0)  T(0)  T(0)  T(0);
        a21   γT    T(0)  T(0)  T(0);
        a31   a32   γT    T(0)  T(0);
        a41   a42   a43   γT    T(0);
        a51   a52   a53   a54   γT
    ]

    b = @SVector [a51, a52, a53, a54, γT]
    c = @SVector [γc, -T2(0.7), T2(0.8), T2(0.924556761814), T2(1)]

    b_embed = @SVector [b1hat2 - a51, b2hat2 - a52, b3hat2 - a53, b4hat2 - a54, -γT]

    return SDIRKTableau(
        A, b, c, γT, 4;
        b_embed = b_embed, embedded_order = 3,
        is_fsal = false, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :default
    )
end

function Kvaerno3Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    # Get correct z-space coefficients from the old tableau
    tab = Kvaerno3Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK: Stage 1 is explicit (FSAL), A[1,:] = 0
    # Stage 2: tmp = uprev + γ*z₁
    # Stage 3: tmp = uprev + a31*z₁ + a32*z₂
    # Stage 4: tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0);
        tab.a31  tab.a32  γT       T(0);
        tab.a41  tab.a42  tab.a43  γT
    ]

    # Stiffly accurate: b = A[4,:] = [a41, a42, a43, γ]
    b = @SVector [tab.a41, tab.a42, tab.a43, γT]
    # c values: master uses c=γ for stage 2
    c = @SVector [T2(0), γc, tab.c3, T2(1)]
    # btilde from old tableau (= bhat - b)
    b_embed = @SVector [tab.btilde1, tab.btilde2, tab.btilde3, tab.btilde4]

    # Predictions: stage 3 from Hermite, stage 4 from yhat
    α_pred = @SMatrix T2[
        0 0 0 0;
        0 0 0 0;
        tab.α31 tab.α32 0 0;
        tab.a31 tab.a32 tab.γ 0
    ]

    return SDIRKTableau(
        A, b, c, γT, 3;
        b_embed = b_embed, embedded_order = 2,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :hermite,
        α_pred = α_pred
    )
end

function KenCarp3Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = KenCarp3Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK: Stage 1 explicit, stages 2-4 implicit
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0);
        tab.a31  tab.a32  γT       T(0);
        tab.a41  tab.a42  tab.a43  γT
    ]

    b = @SVector [tab.a41, tab.a42, tab.a43, γT]
    c = @SVector [T2(0), 2γc, tab.c3, T2(1)]
    b_embed = @SVector [tab.btilde1, tab.btilde2, tab.btilde3, tab.btilde4]

    A_explicit = @SMatrix [
        T(0)      T(0)      T(0)      T(0);
        tab.ea21  T(0)      T(0)      T(0);
        tab.ea31  tab.ea32  T(0)      T(0);
        tab.ea41  tab.ea42  tab.ea43  T(0)
    ]

    b_explicit = @SVector [tab.eb1, tab.eb2, tab.eb3, tab.eb4]
    c_explicit = @SVector [T2(0), 2γc, tab.c3, T2(1)]

    α_pred = @SMatrix T2[
        0 0 0 0;
        0 0 0 0;
        tab.α31 tab.α32 0 0;
        tab.α41 tab.α42 0 0
    ]

    return SDIRKTableau(
        A, b, c, γT, 3;
        b_embed = b_embed, embedded_order = 2,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :kencarp_additive,
        has_additive_splitting = true,
        A_explicit = A_explicit, b_explicit = b_explicit, c_explicit = c_explicit,
        α_pred = α_pred
    )
end

function Kvaerno4Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = Kvaerno4Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 5-stage: stage 1 explicit, stages 2-5 implicit
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT
    ]

    b = @SVector [tab.a51, tab.a52, tab.a53, tab.a54, γT]
    c = @SVector [T2(0), γc, tab.c3, tab.c4, T2(1)]
    b_embed = @SVector [tab.btilde1, tab.btilde2, tab.btilde3, tab.btilde4, tab.btilde5]

    α_pred = @SMatrix T2[
        0 0 0 0 0;
        0 0 0 0 0;
        tab.α31 tab.α32 0 0 0;
        tab.α41 tab.α42 0 0 0;
        tab.a41 tab.a42 tab.a43 tab.γ 0
    ]

    return SDIRKTableau(
        A, b, c, γT, 4;
        b_embed = b_embed, embedded_order = 3,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :hermite,
        α_pred = α_pred
    )
end

function Kvaerno5Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = Kvaerno5Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 7-stage: stage 1 explicit, stages 2-7 implicit
    # Note: a62=0, a72=0 (omitted from old struct)
    A = @SMatrix [
        T(0)     T(0)       T(0)     T(0)     T(0)     T(0)     T(0);
        γT       γT         T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32    γT       T(0)     T(0)     T(0)     T(0);
        tab.a41  tab.a42    tab.a43  γT       T(0)     T(0)     T(0);
        tab.a51  tab.a52    tab.a53  tab.a54  γT       T(0)     T(0);
        tab.a61  T(0)       tab.a63  tab.a64  tab.a65  γT       T(0);
        tab.a71  T(0)       tab.a73  tab.a74  tab.a75  tab.a76  γT
    ]

    b = @SVector [tab.a71, T(0), tab.a73, tab.a74, tab.a75, tab.a76, γT]
    c = @SVector [T2(0), γc, T2(tab.c3), T2(tab.c4), T2(tab.c5), T2(tab.c6), T2(1)]
    b_embed = @SVector [tab.btilde1, T(0), tab.btilde3, tab.btilde4, tab.btilde5, tab.btilde6, tab.btilde7]

    α_pred = @SMatrix T2[
        0 0 0 0 0 0 0;
        0 0 0 0 0 0 0;
        tab.α31 tab.α32 0 0 0 0 0;
        tab.α41 tab.α42 tab.α43 0 0 0 0;
        tab.α51 tab.α52 tab.α53 0 0 0 0;
        tab.α61 tab.α62 tab.α63 0 0 0 0;
        tab.a61 0 tab.a63 tab.a64 tab.a65 tab.γ 0
    ]

    return SDIRKTableau(
        A, b, c, γT, 5;
        b_embed = b_embed, embedded_order = 4,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :hermite,
        α_pred = α_pred
    )
end

function KenCarp4Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = KenCarp4Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 6-stage: stage 1 explicit, stages 2-6 implicit
    # Note: a62=0 (omitted from old struct), btilde2=0 (omitted)
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0);
        tab.a61  T(0)    tab.a63  tab.a64  tab.a65  γT
    ]

    b = @SVector [tab.a61, T(0), tab.a63, tab.a64, tab.a65, γT]
    c = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, T2(1)]
    b_embed = @SVector [tab.btilde1, T(0), tab.btilde3, tab.btilde4, tab.btilde5, tab.btilde6]

    A_explicit = @SMatrix [
        T(0)      T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea21  T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea31  tab.ea32  T(0)      T(0)      T(0)      T(0);
        tab.ea41  tab.ea42  tab.ea43  T(0)      T(0)      T(0);
        tab.ea51  tab.ea52  tab.ea53  tab.ea54  T(0)      T(0);
        tab.ea61  tab.ea62  tab.ea63  tab.ea64  tab.ea65  T(0)
    ]

    b_explicit = @SVector [tab.eb1, T(0), tab.eb3, tab.eb4, tab.eb5, tab.eb6]
    c_explicit = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, T2(1)]

    α_pred = @SMatrix T2[
        0 0 0 0 0 0;
        tab.α21 0 0 0 0 0;
        tab.α31 tab.α32 0 0 0 0;
        tab.α41 tab.α42 0 0 0 0;
        tab.α51 tab.α52 tab.α53 tab.α54 0 0;
        tab.α61 tab.α62 tab.α63 tab.α64 tab.α65 0
    ]

    return SDIRKTableau(
        A, b, c, γT, 4;
        b_embed = b_embed, embedded_order = 3,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :kencarp_additive,
        has_additive_splitting = true,
        A_explicit = A_explicit, b_explicit = b_explicit, c_explicit = c_explicit,
        α_pred = α_pred
    )
end

function KenCarp47Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = KenCarp47Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 7-stage: stage 1 explicit, a71=0, a72=0
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0)     T(0);
        tab.a61  tab.a62  tab.a63  tab.a64  tab.a65  γT       T(0);
        T(0)     T(0)     tab.a73  tab.a74  tab.a75  tab.a76  γT
    ]

    b = @SVector [T(0), T(0), tab.a73, tab.a74, tab.a75, tab.a76, γT]
    c = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, tab.c6, T2(1)]
    b_embed = @SVector [T(0), T(0), tab.btilde3, tab.btilde4, tab.btilde5, tab.btilde6, tab.btilde7]

    A_explicit = @SMatrix [
        T(0)      T(0)      T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea21  T(0)      T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea31  tab.ea32  T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea41  tab.ea42  tab.ea43  T(0)      T(0)      T(0)      T(0);
        tab.ea51  tab.ea52  tab.ea53  tab.ea54  T(0)      T(0)      T(0);
        tab.ea61  tab.ea62  tab.ea63  tab.ea64  tab.ea65  T(0)      T(0);
        tab.ea71  tab.ea72  tab.ea73  tab.ea74  tab.ea75  tab.ea76  T(0)
    ]

    b_explicit = @SVector [T(0), T(0), tab.eb3, tab.eb4, tab.eb5, tab.eb6, tab.eb7]
    c_explicit = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, tab.c6, T2(1)]

    α_pred = @SMatrix T2[
        0 0 0 0 0 0 0;
        tab.α21 0 0 0 0 0 0;
        tab.α31 tab.α32 0 0 0 0 0;
        tab.α41 tab.α42 tab.α43 0 0 0 0;
        tab.α51 tab.α52 0 0 0 0 0;
        tab.α61 tab.α62 tab.α63 0 0 0 0;
        tab.α71 tab.α72 tab.α73 tab.α74 tab.α75 tab.α76 0
    ]

    return SDIRKTableau(
        A, b, c, γT, 4;
        b_embed = b_embed, embedded_order = 3,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :kencarp_additive,
        has_additive_splitting = true,
        A_explicit = A_explicit, b_explicit = b_explicit, c_explicit = c_explicit,
        α_pred = α_pred
    )
end

function KenCarp5Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = KenCarp5Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 8-stage: stage 1 explicit, stages 2-8 implicit
    # Note: a42=0, a82=0, a83=0 (omitted from old struct)
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a41  T(0)    tab.a43  γT       T(0)     T(0)     T(0)     T(0);
        tab.a51  T(0)    tab.a53  tab.a54  γT       T(0)     T(0)     T(0);
        tab.a61  T(0)    tab.a63  tab.a64  tab.a65  γT       T(0)     T(0);
        tab.a71  T(0)    tab.a73  tab.a74  tab.a75  tab.a76  γT       T(0);
        tab.a81  T(0)    T(0)    tab.a84  tab.a85  tab.a86  tab.a87  γT
    ]

    b = @SVector [tab.a81, T(0), T(0), tab.a84, tab.a85, tab.a86, tab.a87, γT]
    c = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7, T2(1)]
    b_embed = @SVector [tab.btilde1, T(0), T(0), tab.btilde4, tab.btilde5, tab.btilde6, tab.btilde7, tab.btilde8]

    A_explicit = @SMatrix [
        T(0)      T(0)      T(0)      T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea21  T(0)      T(0)      T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea31  tab.ea32  T(0)      T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea41  T(0)      tab.ea43  T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea51  T(0)      tab.ea53  tab.ea54  T(0)      T(0)      T(0)      T(0);
        tab.ea61  T(0)      tab.ea63  tab.ea64  tab.ea65  T(0)      T(0)      T(0);
        tab.ea71  T(0)      tab.ea73  tab.ea74  tab.ea75  tab.ea76  T(0)      T(0);
        tab.ea81  T(0)      tab.ea83  tab.ea84  tab.ea85  tab.ea86  tab.ea87  T(0)
    ]

    b_explicit = @SVector [tab.eb1, T(0), T(0), tab.eb4, tab.eb5, tab.eb6, tab.eb7, tab.eb8]
    c_explicit = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7, T2(1)]

    α_pred = @SMatrix T2[
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        tab.α31 tab.α32 0 0 0 0 0 0;
        tab.α41 tab.α42 0 0 0 0 0 0;
        tab.α51 tab.α52 0 0 0 0 0 0;
        tab.α61 tab.α62 0 0 0 0 0 0;
        tab.α71 tab.α72 tab.α73 tab.α74 tab.α75 0 0 0;
        tab.α81 tab.α82 tab.α83 tab.α84 tab.α85 0 0 0
    ]

    return SDIRKTableau(
        A, b, c, γT, 5;
        b_embed = b_embed, embedded_order = 4,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :kencarp_additive,
        has_additive_splitting = true,
        A_explicit = A_explicit, b_explicit = b_explicit, c_explicit = c_explicit,
        α_pred = α_pred
    )
end

function KenCarp58Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = KenCarp58Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 8-stage: stage 1 explicit, a81=0, a82=0
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0)     T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0)     T(0)     T(0);
        tab.a61  tab.a62  tab.a63  tab.a64  tab.a65  γT       T(0)     T(0);
        tab.a71  tab.a72  tab.a73  tab.a74  tab.a75  tab.a76  γT       T(0);
        T(0)     T(0)     tab.a83  tab.a84  tab.a85  tab.a86  tab.a87  γT
    ]

    b = @SVector [T(0), T(0), tab.a83, tab.a84, tab.a85, tab.a86, tab.a87, γT]
    c = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7, T2(1)]
    b_embed = @SVector [T(0), T(0), tab.btilde3, tab.btilde4, tab.btilde5, tab.btilde6, tab.btilde7, tab.btilde8]

    A_explicit = @SMatrix [
        T(0)      T(0)      T(0)      T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea21  T(0)      T(0)      T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea31  tab.ea32  T(0)      T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea41  tab.ea42  tab.ea43  T(0)      T(0)      T(0)      T(0)      T(0);
        tab.ea51  tab.ea52  tab.ea53  tab.ea54  T(0)      T(0)      T(0)      T(0);
        tab.ea61  tab.ea62  tab.ea63  tab.ea64  tab.ea65  T(0)      T(0)      T(0);
        tab.ea71  tab.ea72  tab.ea73  tab.ea74  tab.ea75  tab.ea76  T(0)      T(0);
        tab.ea81  tab.ea82  tab.ea83  tab.ea84  tab.ea85  tab.ea86  tab.ea87  T(0)
    ]

    b_explicit = @SVector [T(0), T(0), tab.eb3, tab.eb4, tab.eb5, tab.eb6, tab.eb7, tab.eb8]
    c_explicit = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7, T2(1)]

    α_pred = @SMatrix T2[
        0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0;
        tab.α31 tab.α32 0 0 0 0 0 0;
        tab.α41 tab.α42 0 0 0 0 0 0;
        tab.α51 tab.α52 0 0 0 0 0 0;
        tab.α61 tab.α62 tab.α63 0 0 0 0 0;
        tab.α71 tab.α72 tab.α73 0 0 0 0 0;
        tab.α81 tab.α82 tab.α83 tab.α84 tab.α85 tab.α86 tab.α87 0
    ]

    return SDIRKTableau(
        A, b, c, γT, 5;
        b_embed = b_embed, embedded_order = 4,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :kencarp_additive,
        has_additive_splitting = true,
        A_explicit = A_explicit, b_explicit = b_explicit, c_explicit = c_explicit,
        α_pred = α_pred
    )
end

function SFSDIRK4Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = SFSDIRK4Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # Pure SDIRK 5-stage (4 implicit + direct solution stage)
    A = @SMatrix [
        γT       T(0)     T(0)     T(0)     T(0);
        tab.a21  γT       T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT
    ]

    b = @SVector [tab.a51, tab.a52, tab.a53, tab.a54, γT]
    c = @SVector [γc, tab.c2, tab.c3, tab.c4, T2(1)]

    return SDIRKTableau(
        A, b, c, γT, 4;
        is_fsal = false, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :default
    )
end

function SFSDIRK5Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = SFSDIRK5Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # Pure SDIRK 6-stage
    A = @SMatrix [
        γT       T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a21  γT       T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0);
        tab.a61  tab.a62  tab.a63  tab.a64  tab.a65  γT
    ]

    b = @SVector [tab.a61, tab.a62, tab.a63, tab.a64, tab.a65, γT]
    c = @SVector [γc, tab.c2, tab.c3, tab.c4, tab.c5, T2(1)]

    return SDIRKTableau(
        A, b, c, γT, 5;
        is_fsal = false, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :default
    )
end

function SFSDIRK6Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = SFSDIRK6Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # Pure SDIRK 7-stage
    A = @SMatrix [
        γT       T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a21  γT       T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0)     T(0);
        tab.a61  tab.a62  tab.a63  tab.a64  tab.a65  γT       T(0);
        tab.a71  tab.a72  tab.a73  tab.a74  tab.a75  tab.a76  γT
    ]

    b = @SVector [tab.a71, tab.a72, tab.a73, tab.a74, tab.a75, tab.a76, γT]
    c = @SVector [γc, tab.c2, tab.c3, tab.c4, tab.c5, tab.c6, T2(1)]

    return SDIRKTableau(
        A, b, c, γT, 6;
        is_fsal = false, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :default
    )
end

function SFSDIRK7Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = SFSDIRK7Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # Pure SDIRK 8-stage
    A = @SMatrix [
        γT       T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a21  γT       T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0)     T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0)     T(0)     T(0);
        tab.a61  tab.a62  tab.a63  tab.a64  tab.a65  γT       T(0)     T(0);
        tab.a71  tab.a72  tab.a73  tab.a74  tab.a75  tab.a76  γT       T(0);
        tab.a81  tab.a82  tab.a83  tab.a84  tab.a85  tab.a86  tab.a87  γT
    ]

    b = @SVector [tab.a81, tab.a82, tab.a83, tab.a84, tab.a85, tab.a86, tab.a87, γT]
    c = @SVector [γc, tab.c2, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7, T2(1)]

    return SDIRKTableau(
        A, b, c, γT, 7;
        is_fsal = false, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :default
    )
end

function SFSDIRK8Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = SFSDIRK8Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # Pure SDIRK 9-stage
    A = @SMatrix [
        γT       T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a21  γT       T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0)     T(0)     T(0)     T(0);
        tab.a61  tab.a62  tab.a63  tab.a64  tab.a65  γT       T(0)     T(0)     T(0);
        tab.a71  tab.a72  tab.a73  tab.a74  tab.a75  tab.a76  γT       T(0)     T(0);
        tab.a81  tab.a82  tab.a83  tab.a84  tab.a85  tab.a86  tab.a87  γT       T(0);
        tab.a91  tab.a92  tab.a93  tab.a94  tab.a95  tab.a96  tab.a97  tab.a98  γT
    ]

    b = @SVector [tab.a91, tab.a92, tab.a93, tab.a94, tab.a95, tab.a96, tab.a97, tab.a98, γT]
    c = @SVector [γc, tab.c2, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7, tab.c8, T2(1)]

    return SDIRKTableau(
        A, b, c, γT, 8;
        is_fsal = false, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :default
    )
end

function ESDIRK54I8L2SATableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = ESDIRK54I8L2SATableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 8-stage: stage 1 explicit
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0)     T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0)     T(0)     T(0);
        tab.a61  tab.a62  tab.a63  tab.a64  tab.a65  γT       T(0)     T(0);
        tab.a71  tab.a72  tab.a73  tab.a74  tab.a75  tab.a76  γT       T(0);
        tab.a81  tab.a82  tab.a83  tab.a84  tab.a85  tab.a86  tab.a87  γT
    ]

    b = @SVector [tab.a81, tab.a82, tab.a83, tab.a84, tab.a85, tab.a86, tab.a87, γT]
    c = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7, T2(1)]
    b_embed = @SVector [tab.btilde1, tab.btilde2, tab.btilde3, tab.btilde4, tab.btilde5, tab.btilde6, tab.btilde7, tab.btilde8]

    return SDIRKTableau(
        A, b, c, γT, 5;
        b_embed = b_embed, embedded_order = 4,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :hermite
    )
end

function ESDIRK436L2SA2Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = ESDIRK436L2SA2Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 6-stage: stage 1 explicit
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0);
        tab.a61  tab.a62  tab.a63  tab.a64  tab.a65  γT
    ]

    b = @SVector [tab.a61, tab.a62, tab.a63, tab.a64, tab.a65, γT]
    c = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, tab.c6]
    b_embed = @SVector [tab.btilde1, tab.btilde2, tab.btilde3, tab.btilde4, tab.btilde5, tab.btilde6]

    return SDIRKTableau(
        A, b, c, γT, 4;
        b_embed = b_embed, embedded_order = 3,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :hermite
    )
end

function ESDIRK437L2SATableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = ESDIRK437L2SATableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 7-stage: stage 1 explicit
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0)     T(0);
        tab.a61  tab.a62  tab.a63  tab.a64  tab.a65  γT       T(0);
        tab.a71  tab.a72  tab.a73  tab.a74  tab.a75  tab.a76  γT
    ]

    b = @SVector [tab.a71, tab.a72, tab.a73, tab.a74, tab.a75, tab.a76, γT]
    c = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7]
    b_embed = @SVector [tab.btilde1, tab.btilde2, tab.btilde3, tab.btilde4, tab.btilde5, tab.btilde6, tab.btilde7]

    return SDIRKTableau(
        A, b, c, γT, 4;
        b_embed = b_embed, embedded_order = 3,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :hermite
    )
end

function ESDIRK547L2SA2Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = ESDIRK547L2SA2Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 7-stage: stage 1 explicit
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0)     T(0);
        tab.a61  tab.a62  tab.a63  tab.a64  tab.a65  γT       T(0);
        tab.a71  tab.a72  tab.a73  tab.a74  tab.a75  tab.a76  γT
    ]

    b = @SVector [tab.a71, tab.a72, tab.a73, tab.a74, tab.a75, tab.a76, γT]
    c = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7]
    b_embed = @SVector [tab.btilde1, tab.btilde2, tab.btilde3, tab.btilde4, tab.btilde5, tab.btilde6, tab.btilde7]

    return SDIRKTableau(
        A, b, c, γT, 5;
        b_embed = b_embed, embedded_order = 4,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :hermite
    )
end

function ESDIRK659L2SATableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = ESDIRK659L2SATableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 9-stage: stage 1 explicit
    # Note: a91=0, a92=0, a93=0 (omitted from old struct)
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0)     T(0)     T(0)     T(0)     T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT       T(0)     T(0)     T(0)     T(0);
        tab.a61  tab.a62  tab.a63  tab.a64  tab.a65  γT       T(0)     T(0)     T(0);
        tab.a71  tab.a72  tab.a73  tab.a74  tab.a75  tab.a76  γT       T(0)     T(0);
        tab.a81  tab.a82  tab.a83  tab.a84  tab.a85  tab.a86  tab.a87  γT       T(0);
        T(0)     T(0)     T(0)     tab.a94  tab.a95  tab.a96  tab.a97  tab.a98  γT
    ]

    b = @SVector [T(0), T(0), T(0), tab.a94, tab.a95, tab.a96, tab.a97, tab.a98, γT]
    c = @SVector [T2(0), 2γc, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7, tab.c8, tab.c9]

    # NOTE: The original btilde values from sdirk_tableaus.jl have sum(btilde) ≈ -1.69,
    # which violates the consistency requirement sum(btilde) = sum(bhat - b) = 0.
    # This appears to be a bug in the original transcription from the Kennedy-Carpenter paper.
    # We correct btilde9 to make sum(btilde) = 0.
    btilde_sum = tab.btilde1 + tab.btilde2 + tab.btilde3 + tab.btilde4 +
        tab.btilde5 + tab.btilde6 + tab.btilde7 + tab.btilde8 + tab.btilde9
    btilde9_corrected = T(tab.btilde9 - btilde_sum)
    b_embed = @SVector [tab.btilde1, tab.btilde2, tab.btilde3, tab.btilde4, tab.btilde5, tab.btilde6, tab.btilde7, tab.btilde8, btilde9_corrected]

    return SDIRKTableau(
        A, b, c, γT, 6;
        b_embed = b_embed, embedded_order = 5,
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :hermite
    )
end

function Hairer4Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = Hairer4Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # Pure SDIRK 5-stage: all stages implicit, A[i,i]=γ
    A = @SMatrix [
        γT       T(0)     T(0)     T(0)     T(0);
        tab.a21  γT       T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT
    ]

    b = @SVector [tab.a51, tab.a52, tab.a53, tab.a54, γT]
    c = @SVector [γc, tab.c2, tab.c3, tab.c4, T2(1)]
    b_embed = @SVector [tab.btilde1, tab.btilde2, tab.btilde3, tab.btilde4, tab.btilde5]

    α_pred = @SMatrix T2[
        0 0 0 0 0;
        tab.α21 0 0 0 0;
        tab.α31 tab.α32 0 0 0;
        tab.α41 0 tab.α43 0 0;
        tab.a41 tab.a42 tab.a43 tab.γ 0
    ]

    return SDIRKTableau(
        A, b, c, γT, 4;
        b_embed = b_embed, embedded_order = 3,
        is_fsal = false, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :hermite,
        α_pred = α_pred
    )
end

function Hairer42Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = Hairer42Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # Pure SDIRK 5-stage
    A = @SMatrix [
        γT       T(0)     T(0)     T(0)     T(0);
        tab.a21  γT       T(0)     T(0)     T(0);
        tab.a31  tab.a32  γT       T(0)     T(0);
        tab.a41  tab.a42  tab.a43  γT       T(0);
        tab.a51  tab.a52  tab.a53  tab.a54  γT
    ]

    b = @SVector [tab.a51, tab.a52, tab.a53, tab.a54, γT]
    c = @SVector [γc, tab.c2, tab.c3, tab.c4, T2(1)]
    b_embed = @SVector [tab.btilde1, tab.btilde2, tab.btilde3, tab.btilde4, tab.btilde5]

    α_pred = @SMatrix T2[
        0 0 0 0 0;
        tab.α21 0 0 0 0;
        tab.α31 tab.α32 0 0 0;
        tab.α41 0 tab.α43 0 0;
        tab.a41 tab.a42 tab.a43 tab.γ 0
    ]

    return SDIRKTableau(
        A, b, c, γT, 4;
        b_embed = b_embed, embedded_order = 3,
        is_fsal = false, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :hermite,
        α_pred = α_pred
    )
end

function CFNLIRK3Tableau_unified(::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = CFNLIRK3Tableau(T, T2)
    γT = T(tab.γ)
    γc = T2(tab.γ)

    # ESDIRK 4-stage: stage 1 explicit
    A = @SMatrix [
        T(0)     T(0)     T(0)     T(0);
        γT       γT       T(0)     T(0);
        tab.a31  tab.a32  γT       T(0);
        tab.a41  tab.a42  tab.a43  γT
    ]

    b = @SVector [tab.a41, tab.a42, tab.a43, γT]
    c = @SVector [T2(0), tab.c2, tab.c3, T2(1)]

    A_explicit = @SMatrix [
        T(0)      T(0)      T(0)      T(0);
        tab.ea21  T(0)      T(0)      T(0);
        tab.ea31  tab.ea32  T(0)      T(0);
        tab.ea41  tab.ea42  tab.ea43  T(0)
    ]

    b_explicit = @SVector [tab.eb1, tab.eb2, tab.eb3, tab.eb4]
    c_explicit = @SVector [T2(0), tab.c2, tab.c3, T2(1)]

    return SDIRKTableau(
        A, b, c, γT, 3;
        is_fsal = true, is_stiffly_accurate = true,
        is_A_stable = true, is_L_stable = true,
        predictor_type = :hermite,
        has_additive_splitting = true,
        A_explicit = A_explicit, b_explicit = b_explicit, c_explicit = c_explicit
    )
end

function get_sdirk_tableau(alg::Symbol, ::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
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
        return Cash4Tableau_unified(T, T2)
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

function get_sdirk_tableau(alg::Cash4, ::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2}
    tab = Cash4Tableau_unified(T, T2)
    if alg.embedding == 2
        b1hat1 = T(1.05646216107052)
        b2hat1 = -T(0.0564621610705236)
        b3hat1 = T(0)
        b4hat1 = T(0)
        b_embed = @SVector [b1hat1 - tab.b[1], b2hat1 - tab.b[2], b3hat1 - tab.b[3], b4hat1 - tab.b[4], -tab.γ]
        return SDIRKTableau(
            tab.A, tab.b, tab.c, tab.γ, 4;
            b_embed = b_embed, embedded_order = 2,
            is_fsal = tab.is_fsal, is_stiffly_accurate = tab.is_stiffly_accurate,
            is_A_stable = tab.is_A_stable, is_L_stable = tab.is_L_stable,
            predictor_type = tab.predictor_type
        )
    end
    return tab
end

get_sdirk_tableau(alg, ::Type{T} = Float64, ::Type{T2} = Float64) where {T, T2} = get_sdirk_tableau(nameof(typeof(alg)), T, T2)
