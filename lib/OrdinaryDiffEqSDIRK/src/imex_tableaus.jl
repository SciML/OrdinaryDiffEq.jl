struct ESDIRKIMEXTableau{T, T2}
    Ai::Matrix{T}
    bi::Vector{T}
    Ae::Matrix{T}
    be::Vector{T}
    c::Vector{T2}
    btilde::Vector{T}
    ebtilde::Union{Vector{T}, Nothing}
    α::Union{Matrix{T2}, Nothing}
    order::Int
    s::Int
end

# Dispatch: each algorithm type maps to its tableau constructor
ESDIRKIMEXTableau(::ARS343, T, T2) = ARS343Tableau(T, T2)
ESDIRKIMEXTableau(::KenCarp3, T, T2) = KenCarp3ESDIRKIMEXTableau(T, T2)
ESDIRKIMEXTableau(::Kvaerno3, T, T2) = Kvaerno3ESDIRKIMEXTableau(T, T2)
ESDIRKIMEXTableau(::Kvaerno4, T, T2) = Kvaerno4ESDIRKIMEXTableau(T, T2)
ESDIRKIMEXTableau(::Kvaerno5, T, T2) = Kvaerno5ESDIRKIMEXTableau(T, T2)
ESDIRKIMEXTableau(::KenCarp4, T, T2) = KenCarp4ESDIRKIMEXTableau(T, T2)
ESDIRKIMEXTableau(::KenCarp5, T, T2) = KenCarp5ESDIRKIMEXTableau(T, T2)
ESDIRKIMEXTableau(::KenCarp47, T, T2) = KenCarp47ESDIRKIMEXTableau(T, T2)
ESDIRKIMEXTableau(::KenCarp58, T, T2) = KenCarp58ESDIRKIMEXTableau(T, T2)

#
# KenCarp3 IMEX Tableau
#

function KenCarp3ESDIRKIMEXTableau(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    γ = convert(T, 0.435866521508459)

    a31 = convert(T, 0.2576482460664272)
    a32 = -convert(T, 0.09351476757488625)
    a41 = convert(T, 0.18764102434672383)
    a42 = -convert(T, 0.595297473576955)
    a43 = convert(T, 0.9717899277217721)

    btilde1 = convert(T, 0.027099261876665316)
    btilde2 = convert(T, 0.11013520969201586)
    btilde3 = -convert(T, 0.10306492520138458)
    btilde4 = -convert(T, 0.0341695463672966)

    c3 = convert(T2, 0.6)
    c2 = 2γ
    θ = c3 / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)
    θ = 1 / c2
    α41 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α42 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)

    ea21 = convert(T, 0.871733043016918)
    ea31 = convert(T, 0.5275890119763004)
    ea32 = convert(T, 0.0724109880236996)
    ea41 = convert(T, 0.3990960076760701)
    ea42 = -convert(T, 0.4375576546135194)
    ea43 = convert(T, 1.0384616469374492)
    eb1 = convert(T, 0.18764102434672383)
    eb2 = -convert(T, 0.595297473576955)
    eb3 = convert(T, 0.9717899277217721)
    eb4 = convert(T, 0.435866521508459)
    ebtilde1 = convert(T, 0.027099261876665316)
    ebtilde2 = convert(T, 0.11013520969201586)
    ebtilde3 = -convert(T, 0.10306492520138458)
    ebtilde4 = -convert(T, 0.0341695463672966)

    s = 4
    Ai = zeros(T, s, s)
    Ai[2, 1] = γ
    Ai[2, 2] = γ
    Ai[3, 1] = a31
    Ai[3, 2] = a32
    Ai[3, 3] = γ
    Ai[4, 1] = a41
    Ai[4, 2] = a42
    Ai[4, 3] = a43
    Ai[4, 4] = γ

    bi_vec = zeros(T, s)
    bi_vec[1] = a41
    bi_vec[2] = a42
    bi_vec[3] = a43
    bi_vec[4] = γ

    Ae = zeros(T, s, s)
    Ae[2, 1] = ea21
    Ae[3, 1] = ea31
    Ae[3, 2] = ea32
    Ae[4, 1] = ea41
    Ae[4, 2] = ea42
    Ae[4, 3] = ea43

    be_vec = zeros(T, s)
    be_vec[1] = eb1
    be_vec[2] = eb2
    be_vec[3] = eb3
    be_vec[4] = eb4

    c_vec = zeros(T2, s)
    c_vec[1] = zero(T2)
    c_vec[2] = convert(T2, 2γ)
    c_vec[3] = c3
    c_vec[4] = one(T2)

    btilde_vec = zeros(T, s)
    btilde_vec[1] = btilde1
    btilde_vec[2] = btilde2
    btilde_vec[3] = btilde3
    btilde_vec[4] = btilde4

    ebtilde_vec = zeros(T, s)
    ebtilde_vec[1] = ebtilde1
    ebtilde_vec[2] = ebtilde2
    ebtilde_vec[3] = ebtilde3
    ebtilde_vec[4] = ebtilde4

    α_mat = zeros(T2, s, s)
    α_mat[3, 1] = α31
    α_mat[3, 2] = α32
    α_mat[4, 1] = α41
    α_mat[4, 2] = α42

    return ESDIRKIMEXTableau(
        Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, ebtilde_vec, α_mat, 3, s
    )
end

function KenCarp3ESDIRKIMEXTableau(T, T2)
    γ = convert(T, 1767732205903 // 4055673282236)

    a31 = convert(T, 2746238789719 // 10658868560708)
    a32 = -convert(T, 640167445237 // 6845629431997)
    a41 = convert(T, 1471266399579 // 7840856788654)
    a42 = -convert(T, 4482444167858 // 7529755066697)
    a43 = convert(T, 11266239266428 // 11593286722821)

    btilde1 = convert(
        T,
        BigInt(681815649026867975666107) //
            BigInt(25159934323302256049469295)
    )
    btilde2 = convert(
        T,
        BigInt(18411887981491912264464127) //
            BigInt(167175311446532472108584143)
    )
    btilde3 = -convert(
        T,
        BigInt(12719313754959329011138489) //
            BigInt(123410692144842870217698057)
    )
    btilde4 = -convert(
        T,
        BigInt(47289384293135913063989) //
            BigInt(1383962894467812063558225)
    )

    c3 = convert(T2, 3 // 5)
    c2 = 2γ
    θ = c3 / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)
    θ = 1 / c2
    α41 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α42 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)

    ea21 = convert(T, 1767732205903 // 2027836641118)
    ea31 = convert(T, 5535828885825 // 10492691773637)
    ea32 = convert(T, 788022342437 // 10882634858940)
    ea41 = convert(T, 6485989280629 // 16251701735622)
    ea42 = -convert(T, 4246266847089 // 9704473918619)
    ea43 = convert(T, 10755448449292 // 10357097424841)
    eb1 = convert(T, 1471266399579 // 7840856788654)
    eb2 = convert(T, -4482444167858 // 7529755066697)
    eb3 = convert(T, 11266239266428 // 11593286722821)
    eb4 = convert(T, 1767732205903 // 4055673282236)
    ebtilde1 = convert(
        T,
        BigInt(681815649026867975666107) //
            BigInt(25159934323302256049469295)
    )
    ebtilde2 = convert(
        T,
        BigInt(18411887981491912264464127) //
            BigInt(167175311446532472108584143)
    )
    ebtilde3 = -convert(
        T,
        BigInt(12719313754959329011138489) //
            BigInt(123410692144842870217698057)
    )
    ebtilde4 = -convert(
        T,
        BigInt(47289384293135913063989) //
            BigInt(1383962894467812063558225)
    )

    s = 4
    Ai = zeros(T, s, s)
    Ai[2, 1] = γ
    Ai[2, 2] = γ
    Ai[3, 1] = a31
    Ai[3, 2] = a32
    Ai[3, 3] = γ
    Ai[4, 1] = a41
    Ai[4, 2] = a42
    Ai[4, 3] = a43
    Ai[4, 4] = γ

    bi_vec = zeros(T, s)
    bi_vec[1] = a41
    bi_vec[2] = a42
    bi_vec[3] = a43
    bi_vec[4] = γ

    Ae = zeros(T, s, s)
    Ae[2, 1] = ea21
    Ae[3, 1] = ea31
    Ae[3, 2] = ea32
    Ae[4, 1] = ea41
    Ae[4, 2] = ea42
    Ae[4, 3] = ea43

    be_vec = zeros(T, s)
    be_vec[1] = eb1
    be_vec[2] = eb2
    be_vec[3] = eb3
    be_vec[4] = eb4

    c_vec = zeros(T2, s)
    c_vec[1] = zero(T2)
    c_vec[2] = convert(T2, 2γ)
    c_vec[3] = c3
    c_vec[4] = one(T2)

    btilde_vec = zeros(T, s)
    btilde_vec[1] = btilde1
    btilde_vec[2] = btilde2
    btilde_vec[3] = btilde3
    btilde_vec[4] = btilde4

    ebtilde_vec = zeros(T, s)
    ebtilde_vec[1] = ebtilde1
    ebtilde_vec[2] = ebtilde2
    ebtilde_vec[3] = ebtilde3
    ebtilde_vec[4] = ebtilde4

    α_mat = zeros(T2, s, s)
    α_mat[3, 1] = α31
    α_mat[3, 2] = α32
    α_mat[4, 1] = α41
    α_mat[4, 2] = α42

    return ESDIRKIMEXTableau(
        Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, ebtilde_vec, α_mat, 3, s
    )
end

#
# ARS343 Tableau
#

function ARS343Tableau(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    γ = convert(T, 0.435866521508459)

    s = 4

    c2 = convert(T2, 0.435866521508459)
    c3 = convert(T2, 0.7179332607542295)
    c4 = one(T2)

    a32_i = convert(T, 0.2820667392457705)

    b3_i = convert(T, -0.644363170684469)
    b2_i = convert(T, 1.20849664917601)

    Ai = zeros(T, s, s)
    Ai[2, 2] = γ
    Ai[3, 2] = a32_i
    Ai[3, 3] = γ
    Ai[4, 2] = b2_i
    Ai[4, 3] = b3_i
    Ai[4, 4] = γ

    bi_vec = T[zero(T), b2_i, b3_i, γ]

    ae21 = convert(T, 0.435866521508459)
    ae31 = convert(T, 0.321278886)
    ae32 = convert(T, 0.3966543748)
    ae41 = -convert(T, 0.105858296)
    ae42 = convert(T, 0.5529291479)
    ae43 = convert(T, 0.5529291479)

    Ae = zeros(T, s, s)
    Ae[2, 1] = ae21
    Ae[3, 1] = ae31
    Ae[3, 2] = ae32
    Ae[4, 1] = ae41
    Ae[4, 2] = ae42
    Ae[4, 3] = ae43

    be_vec = T[zero(T), b2_i, b3_i, γ]

    c_vec = T2[zero(T2), c2, c3, c4]

    btilde_vec = T[zero(T), γ, -2γ, γ]
    ebtilde_vec = T[zero(T), γ, -2γ, γ]

    α_mat = zeros(T2, s, s)

    return ESDIRKIMEXTableau(
        Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, ebtilde_vec, α_mat, 3, s
    )
end

#
# KenCarp4 IMEX Tableau
#

function KenCarp4ESDIRKIMEXTableau(T, T2)
    γ = convert(T2, 1 // 4)

    a31 = convert(T, 8611 // 62500)
    a32 = -convert(T, 1743 // 31250)
    a41 = convert(T, 5012029 // 34652500)
    a42 = -convert(T, 654441 // 2922500)
    a43 = convert(T, 174375 // 388108)
    a51 = convert(T, 15267082809 // 155376265600)
    a52 = -convert(T, 71443401 // 120774400)
    a53 = convert(T, 730878875 // 902184768)
    a54 = convert(T, 2285395 // 8070912)
    a61 = convert(T, 82889 // 524892)
    a63 = convert(T, 15625 // 83664)
    a64 = convert(T, 69875 // 102672)
    a65 = -convert(T, 2260 // 8211)

    btilde1 = convert(T, -31666707 // 9881966720)
    btilde3 = convert(T, 256875 // 105007616)
    btilde4 = convert(T, 2768025 // 128864768)
    btilde5 = -convert(T, 169839 // 3864644)
    btilde6 = convert(T, 5247 // 225920)

    c3 = convert(T2, 83 // 250)
    c4 = convert(T2, 31 // 50)
    c5 = convert(T2, 17 // 20)

    α21 = convert(T2, 2)
    α31 = convert(T2, 42 // 125)
    α32 = convert(T2, 83 // 125)
    α41 = convert(T2, -6 // 25)
    α42 = convert(T2, 31 // 25)
    α51 = convert(T2, 914470432 // 2064665255)
    α52 = convert(T2, 798813 // 724780)
    α53 = convert(T2, -824765625 // 372971788)
    α54 = convert(T2, 49640 // 29791)
    α61 = convert(T2, 288521442795 // 954204491116)
    α62 = convert(T2, 2224881 // 2566456)
    α63 = convert(T2, -1074821875 // 905317354)
    α64 = convert(T2, -3360875 // 8098936)
    α65 = convert(T2, 7040 // 4913)

    ea21 = convert(T, 1 // 2)
    ea31 = convert(T, 13861 // 62500)
    ea32 = convert(T, 6889 // 62500)
    ea41 = -convert(T, 116923316275 // 2393684061468)
    ea42 = -convert(T, 2731218467317 // 15368042101831)
    ea43 = convert(T, 9408046702089 // 11113171139209)
    ea51 = -convert(T, 451086348788 // 2902428689909)
    ea52 = -convert(T, 2682348792572 // 7519795681897)
    ea53 = convert(T, 12662868775082 // 11960479115383)
    ea54 = convert(T, 3355817975965 // 11060851509271)
    ea61 = convert(T, 647845179188 // 3216320057751)
    ea62 = convert(T, 73281519250 // 8382639484533)
    ea63 = convert(T, 552539513391 // 3454668386233)
    ea64 = convert(T, 3354512671639 // 8306763924573)
    ea65 = convert(T, 4040 // 17871)

    eb1 = convert(T, 82889 // 524892)
    eb3 = convert(T, 15625 // 83664)
    eb4 = convert(T, 69875 // 102672)
    eb5 = -convert(T, 2260 // 8211)
    eb6 = convert(T, 1 // 4)

    ebtilde1 = -convert(T, 31666707 // 9881966720)
    ebtilde3 = convert(T, 256875 // 105007616)
    ebtilde4 = convert(T, 2768025 // 128864768)
    ebtilde5 = -convert(T, 169839 // 3864644)
    ebtilde6 = convert(T, 5247 // 225920)

    s = 6
    Ai = zeros(T, s, s)
    Ai[2, 1] = convert(T, γ)
    Ai[2, 2] = convert(T, γ)
    Ai[3, 1] = a31
    Ai[3, 2] = a32
    Ai[3, 3] = convert(T, γ)
    Ai[4, 1] = a41
    Ai[4, 2] = a42
    Ai[4, 3] = a43
    Ai[4, 4] = convert(T, γ)
    Ai[5, 1] = a51
    Ai[5, 2] = a52
    Ai[5, 3] = a53
    Ai[5, 4] = a54
    Ai[5, 5] = convert(T, γ)
    Ai[6, 1] = a61
    Ai[6, 3] = a63
    Ai[6, 4] = a64
    Ai[6, 5] = a65
    Ai[6, 6] = convert(T, γ)

    bi_vec = zeros(T, s)
    bi_vec[1] = a61
    bi_vec[3] = a63
    bi_vec[4] = a64
    bi_vec[5] = a65
    bi_vec[6] = convert(T, γ)

    Ae = zeros(T, s, s)
    Ae[2, 1] = ea21
    Ae[3, 1] = ea31
    Ae[3, 2] = ea32
    Ae[4, 1] = ea41
    Ae[4, 2] = ea42
    Ae[4, 3] = ea43
    Ae[5, 1] = ea51
    Ae[5, 2] = ea52
    Ae[5, 3] = ea53
    Ae[5, 4] = ea54
    Ae[6, 1] = ea61
    Ae[6, 2] = ea62
    Ae[6, 3] = ea63
    Ae[6, 4] = ea64
    Ae[6, 5] = ea65

    be_vec = zeros(T, s)
    be_vec[1] = eb1
    be_vec[3] = eb3
    be_vec[4] = eb4
    be_vec[5] = eb5
    be_vec[6] = eb6

    c_vec = zeros(T2, s)
    c_vec[1] = zero(T2)
    c_vec[2] = convert(T2, 2γ)
    c_vec[3] = c3
    c_vec[4] = c4
    c_vec[5] = c5
    c_vec[6] = one(T2)

    btilde_vec = zeros(T, s)
    btilde_vec[1] = btilde1
    btilde_vec[3] = btilde3
    btilde_vec[4] = btilde4
    btilde_vec[5] = btilde5
    btilde_vec[6] = btilde6

    ebtilde_vec = zeros(T, s)
    ebtilde_vec[1] = ebtilde1
    ebtilde_vec[3] = ebtilde3
    ebtilde_vec[4] = ebtilde4
    ebtilde_vec[5] = ebtilde5
    ebtilde_vec[6] = ebtilde6

    α_mat = zeros(T2, s, s)
    α_mat[2, 1] = α21
    α_mat[3, 1] = α31
    α_mat[3, 2] = α32
    α_mat[4, 1] = α41
    α_mat[4, 2] = α42
    α_mat[5, 1] = α51
    α_mat[5, 2] = α52
    α_mat[5, 3] = α53
    α_mat[5, 4] = α54
    α_mat[6, 1] = α61
    α_mat[6, 2] = α62
    α_mat[6, 3] = α63
    α_mat[6, 4] = α64
    α_mat[6, 5] = α65

    return ESDIRKIMEXTableau(
        Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, ebtilde_vec, α_mat, 4, s
    )
end

#
# KenCarp5 IMEX Tableau
#

function KenCarp5ESDIRKIMEXTableau(T, T2)
    γ = convert(T2, 41 // 200)

    a31 = convert(T, 41 // 400)
    a32 = -convert(T, 567603406766 // 11931857230679)
    a41 = convert(T, 683785636431 // 9252920307686)
    a43 = -convert(T, 110385047103 // 1367015193373)
    a51 = convert(T, 3016520224154 // 10081342136671)
    a53 = convert(T, 30586259806659 // 12414158314087)
    a54 = -convert(T, 22760509404356 // 11113319521817)
    a61 = convert(T, 218866479029 // 1489978393911)
    a63 = convert(T, 638256894668 // 5436446318841)
    a64 = -convert(T, 1179710474555 // 5321154724896)
    a65 = -convert(T, 60928119172 // 8023461067671)
    a71 = convert(T, 1020004230633 // 5715676835656)
    a73 = convert(T, 25762820946817 // 25263940353407)
    a74 = -convert(T, 2161375909145 // 9755907335909)
    a75 = -convert(T, 211217309593 // 5846859502534)
    a76 = -convert(T, 4269925059573 // 7827059040749)
    a81 = -convert(T, 872700587467 // 9133579230613)
    a84 = convert(T, 22348218063261 // 9555858737531)
    a85 = -convert(T, 1143369518992 // 8141816002931)
    a86 = -convert(T, 39379526789629 // 19018526304540)
    a87 = convert(T, 32727382324388 // 42900044865799)

    btilde1 = -convert(T, 360431431567533808054934 // 89473089856732078284381229)
    btilde4 = convert(T, 21220331609936495351431026 // 309921249937726682547321949)
    btilde5 = -convert(T, 42283193605833819490634 // 2144566741190883522084871)
    btilde6 = -convert(T, 21843466548811234473856609 // 296589222149359214696574660)
    btilde7 = convert(T, 3333910710978735057753642 // 199750492790973993533703797)
    btilde8 = convert(T, 45448919757 // 3715198317040)

    c3 = convert(T2, 2935347310677 // 11292855782101)
    c4 = convert(T2, 1426016391358 // 7196633302097)
    c5 = convert(T2, 92 // 100)
    c6 = convert(T2, 24 // 100)
    c7 = convert(T2, 3 // 5)

    α31 = convert(T2, 169472355998441 // 463007087066141)
    α32 = convert(T2, 293534731067700 // 463007087066141)
    α41 = convert(T2, 152460326250177 // 295061965385977)
    α42 = convert(T2, 142601639135800 // 295061965385977)
    α51 = convert(T2, -51 // 41)
    α52 = convert(T2, 92 // 41)
    α61 = convert(T2, 17 // 41)
    α62 = convert(T2, 24 // 41)
    α71 = convert(T2, 13488091065527792 // 122659689776876057)
    α72 = convert(T2, -3214953045 // 3673655312)
    α73 = convert(T2, 550552676519862000 // 151043064207496529)
    α74 = convert(T2, -409689169278408000 // 135215758621947439)
    α75 = convert(T2, 3345 // 12167)
    α81 = convert(T2, 1490668709762032 // 122659689776876057)
    α82 = convert(T2, 5358255075 // 14694621248)
    α83 = convert(T2, -229396948549942500 // 151043064207496529)
    α84 = convert(T2, 170703820532670000 // 135215758621947439)
    α85 = convert(T2, 30275 // 24334)

    ea21 = convert(T, 41 // 100)
    ea31 = convert(T, 367902744464 // 2072280473677)
    ea32 = convert(T, 677623207551 // 8224143866563)
    ea41 = convert(T, 1268023523408 // 10340822734521)
    ea43 = convert(T, 1029933939417 // 13636558850479)
    ea51 = convert(T, 14463281900351 // 6315353703477)
    ea53 = convert(T, 66114435211212 // 5879490589093)
    ea54 = -convert(T, 54053170152839 // 4284798021562)
    ea61 = convert(T, 14090043504691 // 34967701212078)
    ea63 = convert(T, 15191511035443 // 11219624916014)
    ea64 = -convert(T, 18461159152457 // 12425892160975)
    ea65 = -convert(T, 281667163811 // 9011619295870)
    ea71 = convert(T, 19230459214898 // 13134317526959)
    ea73 = convert(T, 21275331358303 // 2942455364971)
    ea74 = -convert(T, 38145345988419 // 4862620318723)
    ea75 = -convert(T, 1 // 8)
    ea76 = -convert(T, 1 // 8)
    ea81 = -convert(T, 19977161125411 // 11928030595625)
    ea83 = -convert(T, 40795976796054 // 6384907823539)
    ea84 = convert(T, 177454434618887 // 12078138498510)
    ea85 = convert(T, 782672205425 // 8267701900261)
    ea86 = -convert(T, 69563011059811 // 9646580694205)
    ea87 = convert(T, 7356628210526 // 4942186776405)

    eb1 = -convert(T, 872700587467 // 9133579230613)
    eb4 = convert(T, 22348218063261 // 9555858737531)
    eb5 = -convert(T, 1143369518992 // 8141816002931)
    eb6 = -convert(T, 39379526789629 // 19018526304540)
    eb7 = convert(T, 32727382324388 // 42900044865799)
    eb8 = convert(T, 41 // 200)

    ebtilde1 = -convert(T, 360431431567533808054934 // 89473089856732078284381229)
    ebtilde4 = convert(T, 21220331609936495351431026 // 309921249937726682547321949)
    ebtilde5 = -convert(T, 42283193605833819490634 // 2144566741190883522084871)
    ebtilde6 = -convert(T, 21843466548811234473856609 // 296589222149359214696574660)
    ebtilde7 = convert(T, 3333910710978735057753642 // 199750492790973993533703797)
    ebtilde8 = convert(T, 45448919757 // 3715198317040)

    s = 8
    Ai = zeros(T, s, s)
    Ai[2, 1] = convert(T, γ)
    Ai[2, 2] = convert(T, γ)
    Ai[3, 1] = a31
    Ai[3, 2] = a32
    Ai[3, 3] = convert(T, γ)
    Ai[4, 1] = a41
    Ai[4, 3] = a43
    Ai[4, 4] = convert(T, γ)
    Ai[5, 1] = a51
    Ai[5, 3] = a53
    Ai[5, 4] = a54
    Ai[5, 5] = convert(T, γ)
    Ai[6, 1] = a61
    Ai[6, 3] = a63
    Ai[6, 4] = a64
    Ai[6, 5] = a65
    Ai[6, 6] = convert(T, γ)
    Ai[7, 1] = a71
    Ai[7, 3] = a73
    Ai[7, 4] = a74
    Ai[7, 5] = a75
    Ai[7, 6] = a76
    Ai[7, 7] = convert(T, γ)
    Ai[8, 1] = a81
    Ai[8, 4] = a84
    Ai[8, 5] = a85
    Ai[8, 6] = a86
    Ai[8, 7] = a87
    Ai[8, 8] = convert(T, γ)

    bi_vec = zeros(T, s)
    bi_vec[1] = a81
    bi_vec[4] = a84
    bi_vec[5] = a85
    bi_vec[6] = a86
    bi_vec[7] = a87
    bi_vec[8] = convert(T, γ)

    Ae = zeros(T, s, s)
    Ae[2, 1] = ea21
    Ae[3, 1] = ea31
    Ae[3, 2] = ea32
    Ae[4, 1] = ea41
    Ae[4, 3] = ea43
    Ae[5, 1] = ea51
    Ae[5, 3] = ea53
    Ae[5, 4] = ea54
    Ae[6, 1] = ea61
    Ae[6, 3] = ea63
    Ae[6, 4] = ea64
    Ae[6, 5] = ea65
    Ae[7, 1] = ea71
    Ae[7, 3] = ea73
    Ae[7, 4] = ea74
    Ae[7, 5] = ea75
    Ae[7, 6] = ea76
    Ae[8, 1] = ea81
    Ae[8, 3] = ea83
    Ae[8, 4] = ea84
    Ae[8, 5] = ea85
    Ae[8, 6] = ea86
    Ae[8, 7] = ea87

    be_vec = zeros(T, s)
    be_vec[1] = eb1
    be_vec[4] = eb4
    be_vec[5] = eb5
    be_vec[6] = eb6
    be_vec[7] = eb7
    be_vec[8] = eb8

    c_vec = zeros(T2, s)
    c_vec[1] = zero(T2)
    c_vec[2] = convert(T2, 2γ)
    c_vec[3] = c3
    c_vec[4] = c4
    c_vec[5] = c5
    c_vec[6] = c6
    c_vec[7] = c7
    c_vec[8] = one(T2)

    btilde_vec = zeros(T, s)
    btilde_vec[1] = btilde1
    btilde_vec[4] = btilde4
    btilde_vec[5] = btilde5
    btilde_vec[6] = btilde6
    btilde_vec[7] = btilde7
    btilde_vec[8] = btilde8

    ebtilde_vec = zeros(T, s)
    ebtilde_vec[1] = ebtilde1
    ebtilde_vec[4] = ebtilde4
    ebtilde_vec[5] = ebtilde5
    ebtilde_vec[6] = ebtilde6
    ebtilde_vec[7] = ebtilde7
    ebtilde_vec[8] = ebtilde8

    α_mat = zeros(T2, s, s)
    α_mat[3, 1] = α31
    α_mat[3, 2] = α32
    α_mat[4, 1] = α41
    α_mat[4, 2] = α42
    α_mat[5, 1] = α51
    α_mat[5, 2] = α52
    α_mat[6, 1] = α61
    α_mat[6, 2] = α62
    α_mat[7, 1] = α71
    α_mat[7, 2] = α72
    α_mat[7, 3] = α73
    α_mat[7, 4] = α74
    α_mat[7, 5] = α75
    α_mat[8, 1] = α81
    α_mat[8, 2] = α82
    α_mat[8, 3] = α83
    α_mat[8, 4] = α84
    α_mat[8, 5] = α85

    return ESDIRKIMEXTableau(
        Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, ebtilde_vec, α_mat, 5, s
    )
end

function ARS343Tableau(T, T2)
    γ = convert(T, 4358665215084590 // 10000000000000000)

    s = 4

    c2 = convert(T2, γ)
    c3 = (one(T2) + convert(T2, γ)) / 2
    c4 = one(T2)

    a32_i = (one(T) - γ) / 2

    b3_i = (one(T) / 2 - 2γ + γ^2) / ((one(T) - γ) / 2)
    b2_i = one(T) - γ - b3_i

    Ai = zeros(T, s, s)
    Ai[2, 2] = convert(T, γ)
    Ai[3, 2] = convert(T, a32_i)
    Ai[3, 3] = convert(T, γ)
    Ai[4, 2] = convert(T, b2_i)
    Ai[4, 3] = convert(T, b3_i)
    Ai[4, 4] = convert(T, γ)

    bi_vec = T[zero(T), convert(T, b2_i), convert(T, b3_i), convert(T, γ)]

    ae21 = convert(T, γ)
    ae31 = convert(T, 3212788860 // 10000000000)
    ae32 = convert(T, 3966543748 // 10000000000)
    ae41 = -convert(T, 1058582960 // 10000000000)
    ae42 = convert(T, 5529291479 // 10000000000)
    ae43 = convert(T, 5529291479 // 10000000000)

    Ae = zeros(T, s, s)
    Ae[2, 1] = ae21
    Ae[3, 1] = ae31
    Ae[3, 2] = ae32
    Ae[4, 1] = ae41
    Ae[4, 2] = ae42
    Ae[4, 3] = ae43

    be_vec = T[zero(T), convert(T, b2_i), convert(T, b3_i), convert(T, γ)]

    c_vec = T2[zero(T2), convert(T2, c2), convert(T2, c3), convert(T2, c4)]

    btilde_vec = T[zero(T), γ, -2γ, γ]
    ebtilde_vec = T[zero(T), γ, -2γ, γ]

    α_mat = zeros(T2, s, s)

    return ESDIRKIMEXTableau(
        Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, ebtilde_vec, α_mat, 3, s
    )
end

#
# Kvaerno3 IMEX Tableau
#

function Kvaerno3ESDIRKIMEXTableau(T, T2)
    γ = convert(T, 0.4358665215)

    a31 = convert(T, 0.490563388419108)
    a32 = convert(T, 0.073570090080892)
    a41 = convert(T, 0.308809969973036)
    a42 = convert(T, 1.490563388254106)
    a43 = -convert(T, 1.235239879727145)

    btilde1 = convert(T, 0.181753418446072)
    btilde2 = convert(T, -1.416993298173214)
    btilde3 = convert(T, 1.671106401227145)
    btilde4 = -convert(T, 0.4358665215)

    c3 = convert(T2, 1)
    c2 = convert(T2, 2) * convert(T2, 0.4358665215)
    θ = c3 / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * convert(T2, 0.4358665215))
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * convert(T2, 0.4358665215))
    α41 = convert(T2, 0.0)
    α42 = convert(T2, 0.0)

    s = 4
    Ai = zeros(T, s, s)
    Ai[2, 1] = γ
    Ai[2, 2] = γ
    Ai[3, 1] = a31
    Ai[3, 2] = a32
    Ai[3, 3] = γ
    Ai[4, 1] = a41
    Ai[4, 2] = a42
    Ai[4, 3] = a43
    Ai[4, 4] = γ

    bi_vec = zeros(T, s)
    bi_vec[1] = a41
    bi_vec[2] = a42
    bi_vec[3] = a43
    bi_vec[4] = γ

    Ae = zeros(T, s, s)
    be_vec = zeros(T, s)

    c_vec = zeros(T2, s)
    c_vec[1] = zero(T2)
    c_vec[2] = c2
    c_vec[3] = c3
    c_vec[4] = one(T2)

    btilde_vec = zeros(T, s)
    btilde_vec[1] = btilde1
    btilde_vec[2] = btilde2
    btilde_vec[3] = btilde3
    btilde_vec[4] = btilde4

    α_mat = zeros(T2, s, s)
    α_mat[3, 1] = α31
    α_mat[3, 2] = α32
    α_mat[4, 1] = α41
    α_mat[4, 2] = α42

    return ESDIRKIMEXTableau(
        Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, nothing, α_mat, 3, s
    )
end

#
# Kvaerno4 IMEX Tableau
#

function Kvaerno4ESDIRKIMEXTableau(T, T2)
    γ = convert(T, 0.4358665215)

    a31 = convert(T, 0.140737774731968)
    a32 = convert(T, -0.108365551378832)
    a41 = convert(T, 0.102399400616089)
    a42 = convert(T, -0.376878452267324)
    a43 = convert(T, 0.838612530151233)
    a51 = convert(T, 0.157024897860995)
    a52 = convert(T, 0.117330441357768)
    a53 = convert(T, 0.61667803039168)
    a54 = convert(T, -0.326899891110444)

    btilde1 = convert(T, -0.054625497244906)
    btilde2 = convert(T, -0.494208893625092)
    btilde3 = convert(T, 0.221934499759553)
    btilde4 = convert(T, 0.762766412610444)
    btilde5 = -convert(T, 0.4358665215)

    c3 = convert(T2, 0.468238744853136)
    c4 = convert(T2, 1)
    c2 = convert(T2, 2) * convert(T2, 0.4358665215)

    α21 = convert(T2, 2)
    α31 = convert(T2, 0.462864521870446)
    α32 = convert(T2, 0.537135478129554)
    α41 = convert(T2, -0.14714018016178376)
    α42 = convert(T2, 1.1471401801617838)

    s = 5
    Ai = zeros(T, s, s)
    Ai[2, 1] = γ
    Ai[2, 2] = γ
    Ai[3, 1] = a31
    Ai[3, 2] = a32
    Ai[3, 3] = γ
    Ai[4, 1] = a41
    Ai[4, 2] = a42
    Ai[4, 3] = a43
    Ai[4, 4] = γ
    Ai[5, 1] = a51
    Ai[5, 2] = a52
    Ai[5, 3] = a53
    Ai[5, 4] = a54
    Ai[5, 5] = γ

    bi_vec = zeros(T, s)
    bi_vec[1] = a51
    bi_vec[2] = a52
    bi_vec[3] = a53
    bi_vec[4] = a54
    bi_vec[5] = γ

    Ae = zeros(T, s, s)
    be_vec = zeros(T, s)

    c_vec = zeros(T2, s)
    c_vec[1] = zero(T2)
    c_vec[2] = c2
    c_vec[3] = c3
    c_vec[4] = c4
    c_vec[5] = one(T2)

    btilde_vec = zeros(T, s)
    btilde_vec[1] = btilde1
    btilde_vec[2] = btilde2
    btilde_vec[3] = btilde3
    btilde_vec[4] = btilde4
    btilde_vec[5] = btilde5

    α_mat = zeros(T2, s, s)
    α_mat[2, 1] = α21
    α_mat[3, 1] = α31
    α_mat[3, 2] = α32
    α_mat[4, 1] = α41
    α_mat[4, 2] = α42

    return ESDIRKIMEXTableau(
        Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, nothing, α_mat, 4, s
    )
end

#
# Kvaerno5 IMEX Tableau
#

function Kvaerno5ESDIRKIMEXTableau(T, T2)
    γ = convert(T, 0.26)

    a31 = convert(T, 0.13)
    a32 = convert(T, 0.84033320996790809)
    a41 = convert(T, 0.22371961478320505)
    a42 = convert(T, 0.47675532319799699)
    a43 = -convert(T, 0.06470895363112615)
    a51 = convert(T, 0.16648564323248321)
    a52 = convert(T, 0.1045001884159172)
    a53 = convert(T, 0.03631482272098715)
    a54 = -convert(T, 0.13090704451073998)
    a61 = convert(T, 0.13855640231268224)
    a63 = -convert(T, 0.04245337201752043)
    a64 = convert(T, 0.02446657898003141)
    a65 = convert(T, 0.61943039072480676)
    a71 = convert(T, 0.13659751177640291)
    a73 = -convert(T, 0.05496908796538376)
    a74 = -convert(T, 0.04118626728321046)
    a75 = convert(T, 0.62993304899016403)
    a76 = convert(T, 0.06962479448202728)

    btilde1 = convert(T, 0.00195889053627933)
    btilde3 = convert(T, 0.01251571594786333)
    btilde4 = convert(T, 0.06565284626324187)
    btilde5 = -convert(T, 0.01050265826535727)
    btilde6 = convert(T, 0.19037520551797272)
    btilde7 = -convert(T, 0.26)

    c3 = convert(T2, 1.230333209967908)
    c4 = convert(T2, 0.895765984350076)
    c5 = convert(T2, 0.436393609858648)
    c6 = convert(T2, 1)
    c2 = convert(T2, 2) * convert(T2, 0.26)

    α21 = convert(T2, 2)
    α31 = convert(T2, -1.366025403784441)
    α32 = convert(T2, 2.3660254037844357)
    α41 = convert(T2, -0.19650552613122207)
    α42 = convert(T2, 0.8113579546496623)
    α43 = convert(T2, 0.38514757148155954)
    α51 = convert(T2, 0.10375304369958693)
    α52 = convert(T2, 0.937994698066431)
    α53 = convert(T2, -0.04174774176601781)
    α61 = convert(T2, -0.17281112873898072)
    α62 = convert(T2, 0.6235784481025847)
    α63 = convert(T2, 0.5492326806363959)

    s = 7
    Ai = zeros(T, s, s)
    Ai[2, 1] = γ
    Ai[2, 2] = γ
    Ai[3, 1] = a31
    Ai[3, 2] = a32
    Ai[3, 3] = γ
    Ai[4, 1] = a41
    Ai[4, 2] = a42
    Ai[4, 3] = a43
    Ai[4, 4] = γ
    Ai[5, 1] = a51
    Ai[5, 2] = a52
    Ai[5, 3] = a53
    Ai[5, 4] = a54
    Ai[5, 5] = γ
    Ai[6, 1] = a61
    Ai[6, 3] = a63
    Ai[6, 4] = a64
    Ai[6, 5] = a65
    Ai[6, 6] = γ
    Ai[7, 1] = a71
    Ai[7, 3] = a73
    Ai[7, 4] = a74
    Ai[7, 5] = a75
    Ai[7, 6] = a76
    Ai[7, 7] = γ

    bi_vec = zeros(T, s)
    bi_vec[1] = a71
    bi_vec[3] = a73
    bi_vec[4] = a74
    bi_vec[5] = a75
    bi_vec[6] = a76
    bi_vec[7] = γ

    Ae = zeros(T, s, s)
    be_vec = zeros(T, s)

    c_vec = zeros(T2, s)
    c_vec[1] = zero(T2)
    c_vec[2] = c2
    c_vec[3] = c3
    c_vec[4] = c4
    c_vec[5] = c5
    c_vec[6] = c6
    c_vec[7] = one(T2)

    btilde_vec = zeros(T, s)
    btilde_vec[1] = btilde1
    btilde_vec[3] = btilde3
    btilde_vec[4] = btilde4
    btilde_vec[5] = btilde5
    btilde_vec[6] = btilde6
    btilde_vec[7] = btilde7

    α_mat = zeros(T2, s, s)
    α_mat[2, 1] = α21
    α_mat[3, 1] = α31
    α_mat[3, 2] = α32
    α_mat[4, 1] = α41
    α_mat[4, 2] = α42
    α_mat[4, 3] = α43
    α_mat[5, 1] = α51
    α_mat[5, 2] = α52
    α_mat[5, 3] = α53
    α_mat[6, 1] = α61
    α_mat[6, 2] = α62
    α_mat[6, 3] = α63

    return ESDIRKIMEXTableau(
        Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, nothing, α_mat, 5, s
    )
end

#
# KenCarp47 IMEX Tableau
#

function KenCarp47ESDIRKIMEXTableau(T, T2)
    γ = convert(T2, 1235 // 10000)

    a31 = convert(T, 624185399699 // 4186980696204)
    a32 = a31
    a41 = convert(T, 1258591069120 // 10082082980243)
    a42 = a41
    a43 = -convert(T, 322722984531 // 8455138723562)
    a51 = -convert(T, 436103496990 // 5971407786587)
    a52 = a51
    a53 = -convert(T, 2689175662187 // 11046760208243)
    a54 = convert(T, 4431412449334 // 12995360898505)
    a61 = -convert(T, 2207373168298 // 14430576638973)
    a62 = a61
    a63 = convert(T, 242511121179 // 3358618340039)
    a64 = convert(T, 3145666661981 // 7780404714551)
    a65 = convert(T, 5882073923981 // 14490790706663)
    a73 = convert(T, 9164257142617 // 17756377923965)
    a74 = -convert(T, 10812980402763 // 74029279521829)
    a75 = convert(T, 1335994250573 // 5691609445217)
    a76 = convert(T, 2273837961795 // 8368240463276)

    btilde3 = convert(T, 216367897668138065439709 // 153341716340757627089664345)
    btilde4 = -convert(T, 1719969231640509698414113 // 303097339249411872572263321)
    btilde5 = convert(T, 33321949854538424751892 // 16748125370719759490730723)
    btilde6 = convert(T, 4033362550194444079469 // 1083063207508329376479196)
    btilde7 = -convert(T, 29 // 20000)

    c3 = convert(T2, 4276536705230 // 10142255878289)
    c4 = convert(T2, 67 // 200)
    c5 = convert(T2, 3 // 40)
    c6 = convert(T2, 7 // 10)

    α21 = convert(T2, 2)
    α31 = -convert(T2, 796131459065721 // 1125899906842624)
    α32 = convert(T2, 961015682954173 // 562949953421312)
    α41 = convert(T2, 139710975840363 // 2251799813685248)
    α42 = convert(T2, 389969885861609 // 1125899906842624)
    α43 = convert(T2, 2664298132243335 // 4503599627370496)
    α51 = convert(T2, 6272219723949193 // 9007199254740992)
    α52 = convert(T2, 2734979530791799 // 9007199254740992)
    α61 = convert(T2, 42616678320173 // 140737488355328)
    α62 = -convert(T2, 2617409280098421 // 1125899906842624)
    α63 = convert(T2, 1701187880189829 // 562949953421312)
    α71 = convert(T2, 4978493057967061 // 2251799813685248)
    α72 = convert(T2, 7230365118049293 // 9007199254740992)
    α73 = -convert(T2, 6826045129237249 // 18014398509481984)
    α74 = -convert(T2, 2388848894891525 // 1125899906842624)
    α75 = -convert(T2, 4796744191239075 // 2251799813685248)
    α76 = convert(T2, 2946706549191323 // 1125899906842624)

    ea21 = convert(T, 247 // 1000)
    ea31 = convert(T, 247 // 4000)
    ea32 = convert(T, 2694949928731 // 7487940209513)
    ea41 = convert(T, 464650059369 // 8764239774964)
    ea42 = convert(T, 878889893998 // 2444806327765)
    ea43 = -convert(T, 952945855348 // 12294611323341)
    ea51 = convert(T, 476636172619 // 8159180917465)
    ea52 = -convert(T, 1271469283451 // 7793814740893)
    ea53 = -convert(T, 859560642026 // 4356155882851)
    ea54 = convert(T, 1723805262919 // 4571918432560)
    ea61 = convert(T, 6338158500785 // 11769362343261)
    ea62 = -convert(T, 4970555480458 // 10924838743837)
    ea63 = convert(T, 3326578051521 // 2647936831840)
    ea64 = -convert(T, 880713585975 // 1841400956686)
    ea65 = -convert(T, 1428733748635 // 8843423958496)
    ea71 = convert(T, 760814592956 // 3276306540349)
    ea72 = convert(T, 760814592956 // 3276306540349)
    ea73 = -convert(T, 47223648122716 // 6934462133451)
    ea74 = convert(T, 71187472546993 // 9669769126921)
    ea75 = -convert(T, 13330509492149 // 9695768672337)
    ea76 = convert(T, 11565764226357 // 8513123442827)

    eb3 = convert(T, 9164257142617 // 17756377923965)
    eb4 = -convert(T, 10812980402763 // 74029279521829)
    eb5 = convert(T, 1335994250573 // 5691609445217)
    eb6 = convert(T, 2273837961795 // 8368240463276)
    eb7 = convert(T, 247 // 2000)

    ebtilde3 = convert(T, 216367897668138065439709 // 153341716340757627089664345)
    ebtilde4 = -convert(T, 1719969231640509698414113 // 303097339249411872572263321)
    ebtilde5 = convert(T, 33321949854538424751892 // 16748125370719759490730723)
    ebtilde6 = convert(T, 4033362550194444079469 // 1083063207508329376479196)
    ebtilde7 = -convert(T, 29 // 20000)

    s = 7
    Ai = zeros(T, s, s)
    Ai[2, 1] = convert(T, γ)
    Ai[2, 2] = convert(T, γ)
    Ai[3, 1] = a31
    Ai[3, 2] = a32
    Ai[3, 3] = convert(T, γ)
    Ai[4, 1] = a41
    Ai[4, 2] = a42
    Ai[4, 3] = a43
    Ai[4, 4] = convert(T, γ)
    Ai[5, 1] = a51
    Ai[5, 2] = a52
    Ai[5, 3] = a53
    Ai[5, 4] = a54
    Ai[5, 5] = convert(T, γ)
    Ai[6, 1] = a61
    Ai[6, 2] = a62
    Ai[6, 3] = a63
    Ai[6, 4] = a64
    Ai[6, 5] = a65
    Ai[6, 6] = convert(T, γ)
    Ai[7, 3] = a73
    Ai[7, 4] = a74
    Ai[7, 5] = a75
    Ai[7, 6] = a76
    Ai[7, 7] = convert(T, γ)

    bi_vec = zeros(T, s)
    bi_vec[3] = a73
    bi_vec[4] = a74
    bi_vec[5] = a75
    bi_vec[6] = a76
    bi_vec[7] = convert(T, γ)

    Ae = zeros(T, s, s)
    Ae[2, 1] = ea21
    Ae[3, 1] = ea31
    Ae[3, 2] = ea32
    Ae[4, 1] = ea41
    Ae[4, 2] = ea42
    Ae[4, 3] = ea43
    Ae[5, 1] = ea51
    Ae[5, 2] = ea52
    Ae[5, 3] = ea53
    Ae[5, 4] = ea54
    Ae[6, 1] = ea61
    Ae[6, 2] = ea62
    Ae[6, 3] = ea63
    Ae[6, 4] = ea64
    Ae[6, 5] = ea65
    Ae[7, 1] = ea71
    Ae[7, 2] = ea72
    Ae[7, 3] = ea73
    Ae[7, 4] = ea74
    Ae[7, 5] = ea75
    Ae[7, 6] = ea76

    be_vec = zeros(T, s)
    be_vec[3] = eb3
    be_vec[4] = eb4
    be_vec[5] = eb5
    be_vec[6] = eb6
    be_vec[7] = eb7

    c_vec = zeros(T2, s)
    c_vec[1] = zero(T2)
    c_vec[2] = convert(T2, 2γ)
    c_vec[3] = c3
    c_vec[4] = c4
    c_vec[5] = c5
    c_vec[6] = c6
    c_vec[7] = one(T2)

    btilde_vec = zeros(T, s)
    btilde_vec[3] = btilde3
    btilde_vec[4] = btilde4
    btilde_vec[5] = btilde5
    btilde_vec[6] = btilde6
    btilde_vec[7] = btilde7

    ebtilde_vec = zeros(T, s)
    ebtilde_vec[3] = ebtilde3
    ebtilde_vec[4] = ebtilde4
    ebtilde_vec[5] = ebtilde5
    ebtilde_vec[6] = ebtilde6
    ebtilde_vec[7] = ebtilde7

    α_mat = zeros(T2, s, s)
    α_mat[2, 1] = α21
    α_mat[3, 1] = α31
    α_mat[3, 2] = α32
    α_mat[4, 1] = α41
    α_mat[4, 2] = α42
    α_mat[4, 3] = α43
    α_mat[5, 1] = α51
    α_mat[5, 2] = α52
    α_mat[6, 1] = α61
    α_mat[6, 2] = α62
    α_mat[6, 3] = α63
    α_mat[7, 1] = α71
    α_mat[7, 2] = α72
    α_mat[7, 3] = α73
    α_mat[7, 4] = α74
    α_mat[7, 5] = α75
    α_mat[7, 6] = α76

    return ESDIRKIMEXTableau(
        Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, ebtilde_vec, α_mat, 4, s
    )
end

#
# KenCarp58 IMEX Tableau
#

function KenCarp58ESDIRKIMEXTableau(T, T2)
    γ = convert(T2, 2 // 9)

    a31 = convert(T, 2366667076620 // 8822750406821)
    a32 = a31
    a41 = -convert(T, 257962897183 // 4451812247028)
    a42 = a41
    a43 = convert(T, 128530224461 // 14379561246022)
    a51 = -convert(T, 486229321650 // 11227943450093)
    a52 = a51
    a53 = -convert(T, 225633144460 // 6633558740617)
    a54 = convert(T, 1741320951451 // 6824444397158)
    a61 = convert(T, 621307788657 // 4714163060173)
    a62 = a61
    a63 = -convert(T, 125196015625 // 3866852212004)
    a64 = convert(T, 940440206406 // 7593089888465)
    a65 = convert(T, 961109811699 // 6734810228204)
    a71 = convert(T, 2036305566805 // 6583108094622)
    a72 = a71
    a73 = -convert(T, 3039402635899 // 4450598839912)
    a74 = -convert(T, 1829510709469 // 31102090912115)
    a75 = -convert(T, 286320471013 // 6931253422520)
    a76 = convert(T, 8651533662697 // 9642993110008)
    a83 = convert(T, 3517720773327 // 20256071687669)
    a84 = convert(T, 4569610470461 // 17934693873752)
    a85 = convert(T, 2819471173109 // 11655438449929)
    a86 = convert(T, 3296210113763 // 10722700128969)
    a87 = -convert(T, 1142099968913 // 5710983926999)

    btilde3 = -convert(T, 18652552508630163520943320 // 168134443655105334713783643)
    btilde4 = convert(T, 141161430501477620145807 // 319735394533244397237135736)
    btilde5 = -convert(T, 207757214437709595456056 // 72283007456311581445415925)
    btilde6 = convert(T, 13674542533282477231637762 // 149163814411398370516486131)
    btilde7 = convert(T, 11939168497868428048898551 // 210101209758476969753215083)
    btilde8 = -convert(T, 1815023333875 // 51666766064334)

    c3 = convert(T2, 6456083330201 // 8509243623797)
    c4 = convert(T2, 1632083962415 // 14158861528103)
    c5 = convert(T2, 6365430648612 // 17842476412687)
    c6 = convert(T2, 18 // 25)
    c7 = convert(T2, 191 // 200)

    α31 = -convert(T2, 796131459065721 // 1125899906842624)
    α32 = convert(T2, 961015682954173 // 562949953421312)
    α41 = convert(T2, 3335563016633385 // 4503599627370496)
    α42 = convert(T2, 2336073221474223 // 9007199254740992)
    α51 = convert(T2, 1777088537295433 // 9007199254740992)
    α52 = convert(T2, 7230110717445555 // 9007199254740992)
    α61 = convert(T2, 305461594360167 // 36028797018963968)
    α62 = convert(T2, 3700851199347703 // 36028797018963968)
    α63 = convert(T2, 8005621056314023 // 9007199254740992)
    α71 = convert(T2, 247009276011491 // 9007199254740992)
    α72 = -convert(T2, 6222030107065861 // 9007199254740992)
    α73 = convert(T2, 1872777510724421 // 1125899906842624)
    α81 = convert(T2, 180631849429283 // 36028797018963968)
    α82 = -convert(T2, 3454740038041085 // 36028797018963968)
    α83 = convert(T2, 476708848972457 // 2251799813685248)
    α84 = convert(T2, 5255799236757313 // 288230376151711744)
    α85 = convert(T2, 3690914796734375 // 288230376151711744)
    α86 = -convert(T2, 5010195363762467 // 18014398509481984)
    α87 = convert(T2, 5072201887169367 // 4503599627370496)

    ea21 = convert(T, 4 // 9)
    ea31 = convert(T, 1 // 9)
    ea32 = convert(T, 1183333538310 // 1827251437969)
    ea41 = convert(T, 895379019517 // 9750411845327)
    ea42 = convert(T, 477606656805 // 13473228687314)
    ea43 = -convert(T, 112564739183 // 9373365219272)
    ea51 = -convert(T, 4458043123994 // 13015289567637)
    ea52 = -convert(T, 2500665203865 // 9342069639922)
    ea53 = convert(T, 983347055801 // 8893519644487)
    ea54 = convert(T, 2185051477207 // 2551468980502)
    ea61 = -convert(T, 167316361917 // 17121522574472)
    ea62 = convert(T, 1605541814917 // 7619724128744)
    ea63 = convert(T, 991021770328 // 13052792161721)
    ea64 = convert(T, 2342280609577 // 11279663441611)
    ea65 = convert(T, 3012424348531 // 12792462456678)
    ea71 = convert(T, 6680998715867 // 14310383562358)
    ea72 = convert(T, 5029118570809 // 3897454228471)
    ea73 = convert(T, 2415062538259 // 6382199904604)
    ea74 = -convert(T, 3924368632305 // 6964820224454)
    ea75 = -convert(T, 4331110370267 // 15021686902756)
    ea76 = -convert(T, 3944303808049 // 11994238218192)
    ea81 = convert(T, 2193717860234 // 3570523412979)
    ea82 = convert(T, 2193717860234 // 3570523412979)
    ea83 = convert(T, 5952760925747 // 18750164281544)
    ea84 = -convert(T, 4412967128996 // 6196664114337)
    ea85 = convert(T, 4151782504231 // 36106512998704)
    ea86 = convert(T, 572599549169 // 6265429158920)
    ea87 = -convert(T, 457874356192 // 11306498036315)

    eb3 = convert(T, 3517720773327 // 20256071687669)
    eb4 = convert(T, 4569610470461 // 17934693873752)
    eb5 = convert(T, 2819471173109 // 11655438449929)
    eb6 = convert(T, 3296210113763 // 10722700128969)
    eb7 = -convert(T, 1142099968913 // 5710983926999)
    eb8 = convert(T, 2 // 9)

    ebtilde3 = -convert(T, 18652552508630163520943320 // 168134443655105334713783643)
    ebtilde4 = convert(T, 141161430501477620145807 // 319735394533244397237135736)
    ebtilde5 = -convert(T, 207757214437709595456056 // 72283007456311581445415925)
    ebtilde6 = convert(T, 13674542533282477231637762 // 149163814411398370516486131)
    ebtilde7 = convert(T, 11939168497868428048898551 // 210101209758476969753215083)
    ebtilde8 = -convert(T, 1815023333875 // 51666766064334)

    s = 8
    Ai = zeros(T, s, s)
    Ai[2, 1] = convert(T, γ)
    Ai[2, 2] = convert(T, γ)
    Ai[3, 1] = a31
    Ai[3, 2] = a32
    Ai[3, 3] = convert(T, γ)
    Ai[4, 1] = a41
    Ai[4, 2] = a42
    Ai[4, 3] = a43
    Ai[4, 4] = convert(T, γ)
    Ai[5, 1] = a51
    Ai[5, 2] = a52
    Ai[5, 3] = a53
    Ai[5, 4] = a54
    Ai[5, 5] = convert(T, γ)
    Ai[6, 1] = a61
    Ai[6, 2] = a62
    Ai[6, 3] = a63
    Ai[6, 4] = a64
    Ai[6, 5] = a65
    Ai[6, 6] = convert(T, γ)
    Ai[7, 1] = a71
    Ai[7, 2] = a72
    Ai[7, 3] = a73
    Ai[7, 4] = a74
    Ai[7, 5] = a75
    Ai[7, 6] = a76
    Ai[7, 7] = convert(T, γ)
    Ai[8, 3] = a83
    Ai[8, 4] = a84
    Ai[8, 5] = a85
    Ai[8, 6] = a86
    Ai[8, 7] = a87
    Ai[8, 8] = convert(T, γ)

    bi_vec = zeros(T, s)
    bi_vec[3] = a83
    bi_vec[4] = a84
    bi_vec[5] = a85
    bi_vec[6] = a86
    bi_vec[7] = a87
    bi_vec[8] = convert(T, γ)

    Ae = zeros(T, s, s)
    Ae[2, 1] = ea21
    Ae[3, 1] = ea31
    Ae[3, 2] = ea32
    Ae[4, 1] = ea41
    Ae[4, 2] = ea42
    Ae[4, 3] = ea43
    Ae[5, 1] = ea51
    Ae[5, 2] = ea52
    Ae[5, 3] = ea53
    Ae[5, 4] = ea54
    Ae[6, 1] = ea61
    Ae[6, 2] = ea62
    Ae[6, 3] = ea63
    Ae[6, 4] = ea64
    Ae[6, 5] = ea65
    Ae[7, 1] = ea71
    Ae[7, 2] = ea72
    Ae[7, 3] = ea73
    Ae[7, 4] = ea74
    Ae[7, 5] = ea75
    Ae[7, 6] = ea76
    Ae[8, 1] = ea81
    Ae[8, 2] = ea82
    Ae[8, 3] = ea83
    Ae[8, 4] = ea84
    Ae[8, 5] = ea85
    Ae[8, 6] = ea86
    Ae[8, 7] = ea87

    be_vec = zeros(T, s)
    be_vec[3] = eb3
    be_vec[4] = eb4
    be_vec[5] = eb5
    be_vec[6] = eb6
    be_vec[7] = eb7
    be_vec[8] = eb8

    c_vec = zeros(T2, s)
    c_vec[1] = zero(T2)
    c_vec[2] = convert(T2, 2γ)
    c_vec[3] = c3
    c_vec[4] = c4
    c_vec[5] = c5
    c_vec[6] = c6
    c_vec[7] = c7
    c_vec[8] = one(T2)

    btilde_vec = zeros(T, s)
    btilde_vec[3] = btilde3
    btilde_vec[4] = btilde4
    btilde_vec[5] = btilde5
    btilde_vec[6] = btilde6
    btilde_vec[7] = btilde7
    btilde_vec[8] = btilde8

    ebtilde_vec = zeros(T, s)
    ebtilde_vec[3] = ebtilde3
    ebtilde_vec[4] = ebtilde4
    ebtilde_vec[5] = ebtilde5
    ebtilde_vec[6] = ebtilde6
    ebtilde_vec[7] = ebtilde7
    ebtilde_vec[8] = ebtilde8

    α_mat = zeros(T2, s, s)
    α_mat[3, 1] = α31
    α_mat[3, 2] = α32
    α_mat[4, 1] = α41
    α_mat[4, 2] = α42
    α_mat[5, 1] = α51
    α_mat[5, 2] = α52
    α_mat[6, 1] = α61
    α_mat[6, 2] = α62
    α_mat[6, 3] = α63
    α_mat[7, 1] = α71
    α_mat[7, 2] = α72
    α_mat[7, 3] = α73
    α_mat[8, 1] = α81
    α_mat[8, 2] = α82
    α_mat[8, 3] = α83
    α_mat[8, 4] = α84
    α_mat[8, 5] = α85
    α_mat[8, 6] = α86
    α_mat[8, 7] = α87

    return ESDIRKIMEXTableau(
        Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, ebtilde_vec, α_mat, 5, s
    )
end
