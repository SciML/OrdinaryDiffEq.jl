struct RadauIIA3Tableau{T, T2}
    T11::T
    T12::T
    T21::T
    T22::T
    TI11::T
    TI12::T
    TI21::T
    TI22::T
    c1::T2
    c2::T2
    α::T
    β::T
    e1::T
    e2::T
end

function RadauIIA3Tableau(T, T2)
    T11 = T(0.10540925533894596)
    T12 = T(-0.29814239699997197)
    T21 = T(0.9486832980505138)
    T22 = T(0.0)
    TI11 = T(0.0)
    TI12 = T(1.0540925533894598)
    TI21 = T(-3.3541019662496843)
    TI22 = T(0.3726779962499649)
    c1 = T2(1 / 3)
    c2 = T2(1.0)
    α = T(2.0)
    β = T(-sqrt(2))
    e1 = T(1 / 4)
    e2 = T(-1 / 4)
    RadauIIA3Tableau{T, T2}(T11, T12, T21, T22,
        TI11, TI12, TI21, TI22,
        c1, c2, α, β, e1, e2)
end

struct RadauIIA5Tableau{T, T2}
    T11::T
    T12::T
    T13::T
    T21::T
    T22::T
    T23::T
    T31::T
    #T32::T
    #T33::T
    TI11::T
    TI12::T
    TI13::T
    TI21::T
    TI22::T
    TI23::T
    TI31::T
    TI32::T
    TI33::T
    c1::T2
    c2::T2
    #c3::T2
    γ::T
    α::T
    β::T
    e1::T
    e2::T
    e3::T
end

# inv(T) * inv(A) * T = [γ  0  0
#                        0  α -β
#                        0  β  α]
function RadauIIA5Tableau(T, T2)
    T11 = convert(T, 9.1232394870892942792e-02)
    T12 = convert(T, -0.14125529502095420843e0)
    T13 = convert(T, -3.0029194105147424492e-02)
    T21 = convert(T, 0.24171793270710701896e0)
    T22 = convert(T, 0.20412935229379993199e0)
    T23 = convert(T, 0.38294211275726193779e0)
    T31 = convert(T, 0.96604818261509293619e0)
    TI11 = convert(T, 4.3255798900631553510e0)
    TI12 = convert(T, 0.33919925181580986954e0)
    TI13 = convert(T, 0.54177053993587487119e0)
    TI21 = convert(T, -4.1787185915519047273e0)
    TI22 = convert(T, -0.32768282076106238708e0)
    TI23 = convert(T, 0.47662355450055045196e0)
    TI31 = convert(T, -0.50287263494578687595e0)
    TI32 = convert(T, 2.5719269498556054292e0)
    TI33 = convert(T, -0.59603920482822492497e0)

    sqrt6 = sqrt(6)

    c1 = convert(T2, (4 - sqrt6) / 10)
    c2 = convert(T2, (4 + sqrt6) / 10)

    cbrt9 = cbrt(9)

    γ′ = convert(T, (6.0 + cbrt9 * (cbrt9 - 1)) / 30) # eigval of `A`
    α′ = convert(T, (12.0 - cbrt9 * (cbrt9 - 1)) / 60) # eigval of `A`
    β′ = convert(T, cbrt9 * (cbrt9 + 1) * sqrt(3) / 60) # eigval of `A`
    scale = α′^2 + β′^2
    γ = inv(γ′)  # eigval of `inv(A)`
    α = α′ / scale # eigval of `inv(A)`
    β = β′ / scale # eigval of `inv(A)`

    e1 = convert(T, -(13 + 7 * sqrt6) / 3)
    e2 = convert(T, (-13 + 7 * sqrt6) / 3)
    e3 = convert(T, -1 / 3)
    RadauIIA5Tableau{T, T2}(T11, T12, T13, T21, T22, T23, T31, #= T33 = 0 =#
        TI11, TI12, TI13, TI21, TI22, TI23, TI31, TI32, TI33,
        c1, c2, #= c3 = 1 =#
        γ, α, β,
        e1, e2, e3)
end

struct RadauIIA9Tableau{T, T2}
    T11::T
    T12::T
    T13::T
    T14::T
    T15::T
    T21::T
    T22::T
    T23::T
    T24::T
    T25::T
    T31::T
    T32::T
    T33::T
    T34::T
    T35::T
    T41::T
    T42::T
    T43::T
    T44::T
    T45::T
    T51::T
    #T52::T
    #T53::T
    #T54::T
    #T55::T
    TI11::T
    TI12::T
    TI13::T
    TI14::T
    TI15::T
    TI21::T
    TI22::T
    TI23::T
    TI24::T
    TI25::T
    TI31::T
    TI32::T
    TI33::T
    TI34::T
    TI35::T
    TI41::T
    TI42::T
    TI43::T
    TI44::T
    TI45::T
    TI51::T
    TI52::T
    TI53::T
    TI54::T
    TI55::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    γ::T
    α1::T
    β1::T
    α2::T
    β2::T
    e1::T
    e2::T
    e3::T
    e4::T
    e5::T
end

function RadauIIA9Tableau(T, T2)
    T11 = convert(T, -1.251758622050104589014e-2)
    T12 = convert(T, -1.024204781790882707009e-2)
    T13 = convert(T, 4.767387729029572386318e-2)
    T14 = convert(T, -1.147851525522951470794e-2)
    T15 = convert(T, -1.401985889287541028108e-2)
    T21 = convert(T, -1.491670151895382429004e-3)
    T22 = convert(T, 5.017286451737105816299e-2)
    T23 = convert(T, -9.433181918161143698066e-2)
    T24 = convert(T, -7.668830749180162885157e-3)
    T25 = convert(T, 2.470857842651852681253e-2)
    T31 = convert(T, -7.298187638808714862266e-2)
    T32 = convert(T, -2.305395340434179467214e-1)
    T33 = convert(T, 1.027030453801258997922e-1)
    T34 = convert(T, 1.939846399882895091122e-2)
    T35 = convert(T, 8.180035370375117083639e-2)
    T41 = convert(T, -3.800914400035681041264e-1)
    T42 = convert(T, 3.778939022488612495439e-1)
    T43 = convert(T, 4.667441303324943592896e-1)
    T44 = convert(T, 4.076011712801990666217e-1)
    T45 = convert(T, 1.996824278868025259365e-1)
    T51 = convert(T, -9.219789736812104884883e-1)
    #T52 = convert(T, 1.000000000000000000000e0)
    #T53 = convert(T, 0.0000000000000000000000e0)
    #T54 = convert(T, 1.000000000000000000000e0)
    #T55 = convert(T, 0.0000000000000000000000e0)
    TI11 = convert(T, -3.004156772154440162771e1)
    TI12 = convert(T, -1.386510785627141316518e1)
    TI13 = convert(T, -3.480002774795185561828e0)
    TI14 = convert(T, 1.032008797825263422771e0)
    TI15 = convert(T, -8.043030450739899174753e-1)
    TI21 = convert(T, 5.344186437834911598895e0)
    TI22 = convert(T, 4.593615567759161004454e0)
    TI23 = convert(T, -3.036360323459424298646e0)
    TI24 = convert(T, 1.050660190231458863860e0)
    TI25 = convert(T, -2.727786118642962705386e-1)
    TI31 = convert(T, 3.748059807439804860051e0)
    TI32 = convert(T, -3.984965736343884667252e0)
    TI33 = convert(T, -1.044415641608018792942e0)
    TI34 = convert(T, 1.184098568137948487231e0)
    TI35 = convert(T, -4.499177701567803688988e-1)
    TI41 = convert(T, -3.304188021351900000806e1)
    TI42 = convert(T, -1.737695347906356701945e1)
    TI43 = convert(T, -1.721290632540055611515e-1)
    TI44 = convert(T, -9.916977798254264258817e-2)
    TI45 = convert(T, 5.312281158383066671849e-1)
    TI51 = convert(T, -8.611443979875291977700e0)
    TI52 = convert(T, 9.699991409528808231336e0)
    TI53 = convert(T, 1.914728639696874284851e0)
    TI54 = convert(T, 2.418692006084940026427e0)
    TI55 = convert(T, -1.047463487935337418694e0)

    c1 = convert(T2, 5.710419611451768219312e-2)
    c2 = convert(T2, 2.768430136381238276800e-1)
    c3 = convert(T2, 5.835904323689168200567e-1)
    c4 = convert(T2, 8.602401356562194478479e-1)
    #= c5 = convert(T2, 1) =#

    γ = convert(T, 6.286704751729276645173e0)
    α1 = convert(T, 3.655694325463572258243e0)
    β1 = convert(T, 6.543736899360077294021e0)
    α2 = convert(T, 5.700953298671789419170e0)
    β2 = convert(T, 3.210265600308549888425e0)

    e1 = convert(T, -2.778093394406463730479e1)
    e2 = convert(T, 3.641478498049213152712e0)
    e3 = convert(T, -1.252547721169118720491e0)
    e4 = convert(T, 5.920031671845428725662e-1)
    e5 = convert(T, -2.000000000000000000000e-1)

    RadauIIA9Tableau{T, T2}(T11, T12, T13, T14, T15,
        T21, T22, T23, T24, T25, T31, T32, T33, T34, T35,
        T41, T42, T43, T44, T45, T51, #=T52, T53, T54, T55=#
        TI11, TI12, TI13, TI14, TI15, TI21, TI22, TI23, TI24, TI25,
        TI31, TI32, TI33, TI34, TI35, TI41, TI42, TI43, TI44, TI45,
        TI51, TI52, TI53, TI54, TI55,
        c1, c2, c3, c4, #= c5 = 1 =#
        γ, α1, β1, α2, β2,
        e1, e2, e3, e4, e5)
end

struct RadauIIATableau{T1, T2}
    T::Matrix{T1}
    TI::Matrix{T1}
    c::Vector{T2}
    γ::T1
    α::Vector{T1}
    β::Vector{T1}
    e::Vector{T1}
end

import LinearAlgebra: eigen
import FastGaussQuadrature: gaussradau

function RadauIIATableau{T1, T2}(tab::RadauIIATableau{T1, T2}) where {T1, T2}
    RadauIIATableau{T1, T2}(tab.T, tab.TI, tab.c, tab.γ, tab.α, tab.β, tab.e)
end

function RadauIIATableau(T1, T2, num_stages::Int)
    tab = get(RadauIIATableauCache, (T1, T2, num_stages)) do
        tab = generateRadauTableau(T1, T2, num_stages)
        RadauIIATableauCache[T1, T2, num_stages] = tab
        tab
    end
    return RadauIIATableau{T1, T2}(tab)
end

function generateRadauTableau(T1, T2, num_stages::Int)
    c = reverse!(1 .- gaussradau(num_stages, T1)[1]) ./ 2
    if T1 == T2
        c2 = c
    else
        c2 = reverse!(1 .- gaussradau(num_stages, T2)[1]) ./ 2
    end

    c_powers = Matrix{T1}(undef, num_stages, num_stages)
    for i in 1:num_stages
        c_powers[i, 1] = 1
        for j in 2:num_stages
            c_powers[i, j] = c[i] * c_powers[i, j - 1]
        end
    end
    c_q = Matrix{T1}(undef, num_stages, num_stages)
    for i in 1:num_stages
        for j in 1:num_stages
            c_q[i, j] = c_powers[i, j] * c[i] / j
        end
    end
    a = c_q / c_powers

    local eigval, eigvec
    try
        eigval, eigvec = eigen(a)
    catch
        throw(ArgumentError("Solving ODEs with AdaptiveRadau with $T1 eltype and max_order >=17 requires loading GenericSchur.jl"))
    end
    # α, β, and γ come from eigvals(inv(a)) which are equal to inv.(eivals(a))
    eigval .= inv.(eigval)
    α = [real(eigval[i]) for i in 1:2:(num_stages - 1)]
    β = [imag(eigval[i]) for i in 1:2:(num_stages - 1)]
    γ = real(eigval[num_stages])

    T = Matrix{T1}(undef, num_stages, num_stages)
    @views for i in 2:2:num_stages
        T[:, i] .= real.(eigvec[:, i] ./ eigvec[num_stages, i])
        T[:, i + 1] .= imag.(eigvec[:, i] ./ eigvec[num_stages, i])
    end
    @views T[:, 1] .= real.(eigvec[:, num_stages])
    TI = inv(T)
    # TODO: figure out why all the order conditions are the same
    A = c_powers' ./ (1:num_stages)
    # TODO: figure out why these are the right b
    b = vcat(-(num_stages)^2, -0.5, zeros(num_stages - 2))
    e = A \ b
    tab = RadauIIATableau{T1, T2}(T, TI, c2, γ, α, β, e)
end

const RadauIIATableauCache = Dict{
    Tuple{Type, Type, Int}, RadauIIATableau{T1, T2} where {T1, T2}}(
    (Float64, Float64, 3) => generateRadauTableau(Float64, Float64, 3),
    (Float64, Float64, 5) => generateRadauTableau(Float64, Float64, 5),
    (Float64, Float64, 7) => generateRadauTableau(Float64, Float64, 7))
