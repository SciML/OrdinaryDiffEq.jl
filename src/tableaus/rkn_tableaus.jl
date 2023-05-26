struct Nystrom5ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2
    c7::T2
    a21::T
    a31::T
    a32::T
    a41::T
    #a42::T
    a43::T
    a51::T
    a52::T
    a53::T
    a54::T
    a61::T
    a62::T
    a63::T
    a64::T
    #a65::T
    a71::T
    #a72::T
    a73::T
    a74::T
    a75::T
    #a76::T
    abar21::T
    abar31::T
    abar32::T
    abar41::T
    abar42::T
    abar43::T
    abar51::T
    abar52::T
    abar53::T
    abar54::T
    abar61::T
    abar62::T
    abar63::T
    abar64::T
    abar65::T
    abar71::T
    #abar72::T
    abar73::T
    abar74::T
    abar75::T
    abar76::T
    b1::T
    #b2::T
    b3::T
    b4::T
    b5::T
    #b6::T
    #b7::T
    bbar1::T
    #bbar2::T
    bbar3::T
    bbar4::T
    bbar5::T
    bbar6::T
    #bbar7::T
end

function Nystrom5ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    c1 = convert(T2, 1.0)
    c2 = convert(T2, 0.2051282051282051)
    c3 = convert(T2, 0.307692307692308)
    c4 = convert(T2, 0.8333333333333334)
    c5 = convert(T2, 0.9148936170212766)
    c6 = convert(T2, 1.0)
    c7 = convert(T2, 1.0)
    a21 = convert(T, 0.0210387902695595)
    a31 = convert(T, 0.0236686390532544)
    a32 = convert(T, 0.0236686390532544)
    a41 = convert(T, 0.0337577160493827)
    #a42 = convert(T, 0.0)
    a43 = convert(T, 0.31346450617284)
    a51 = convert(T, -0.060954498687391)
    a52 = convert(T, 0.145802190676171)
    a53 = convert(T, 0.34307079575472)
    a54 = convert(T, -0.0094033225103623)
    a61 = convert(T, -0.135737843227489)
    a62 = convert(T, 0.352968246663599)
    a63 = convert(T, 0.268963798128547)
    a64 = convert(T, 0.0138057984353428)
    #a65 = convert(T, 0.0)
    a71 = convert(T, 0.0933527131782946)
    #a72 = convert(T, 0.0)
    a73 = convert(T, 0.319562323318651)
    a74 = convert(T, 0.138960763520679)
    a75 = convert(T, -0.0518758000176242)
    #a76 = convert(T, 0.0)
    abar21 = convert(T, 0.205128205128205)
    abar31 = convert(T, 0.0769230769230769)
    abar32 = convert(T, 0.230769230769231)
    abar41 = convert(T, 1.06843171296296)
    abar42 = convert(T, -4.09071180555556)
    abar43 = convert(T, 3.85561342592593)
    abar51 = convert(T, 2.44433442449694)
    abar52 = convert(T, -9.53630178503099)
    abar53 = convert(T, 8.17612078274186)
    abar54 = convert(T, -0.169259805186522)
    abar61 = convert(T, 2.9748519087319)
    abar62 = convert(T, -11.3596698113208)
    abar63 = convert(T, 9.45765911709872)
    abar64 = convert(T, 0.162068068588807)
    abar65 = convert(T, -0.234909283098676)
    abar71 = convert(T, 0.0933527131782946)
    #abar72 = convert(T, 0.0)
    abar73 = convert(T, 0.461590022571385)
    abar74 = convert(T, 0.833764581124072)
    abar75 = convert(T, -0.609540650207085)
    abar76 = convert(T, 0.220833333333333)
    b1 = convert(T, 0.0933527131782946)
    #b2 = convert(T, 0.0)
    b3 = convert(T, 0.319562323318651)
    b4 = convert(T, 0.138960763520679)
    b5 = convert(T, -0.0518758000176242)
    #b6 = convert(T, 0.0)
    #b7 = convert(T, 0.0)
    bbar1 = convert(T, 0.0933527131782946)
    bbar2 = convert(T, 0.0)
    bbar3 = convert(T, 0.461590022571385)
    bbar4 = convert(T, 0.833764581124072)
    bbar5 = convert(T, -0.609540650207085)
    bbar6 = convert(T, 0.220833333333333)
    #bbar7 = convert(T, 0.0)
    Nystrom5ConstantCache(c1, c2, c3, c4, c5, c6, c7, a21, a31, a32, a41, a43, a51,
                          a52, a53, a54, a61, a62, a63, a64, a71, a73, a74, a75,
                          abar21, abar31, abar32, abar41, abar42, abar43, abar51,
                          abar52, abar53, abar54, abar61, abar62, abar63, abar64, abar65,
                          abar71, abar73, abar74, abar75, abar76, b1, b3, b4,
                          b5, bbar1, bbar3, bbar4, bbar5, bbar6)
end

function Nystrom5ConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 1 // 1)
    c2 = convert(T2, 8 // 39)
    c3 = convert(T2, 4 // 13)
    c4 = convert(T2, 5 // 6)
    c5 = convert(T2, 43 // 47)
    c6 = convert(T2, 1 // 1) # 36463 // 36464
    c7 = convert(T2, 1 // 1)
    a21 = convert(T, 32 // 1521)
    a31 = convert(T, 4 // 169)
    a32 = convert(T, 4 // 169)
    a41 = convert(T, 175 // 5184)
    #a42 = convert(T, 0 // 1)
    a43 = convert(T, 1625 // 5184)
    a51 = convert(T, -342497279 // 5618900760)
    a52 = convert(T, 6827067 // 46824173)
    a53 = convert(T, 35048741 // 102161832)
    a54 = convert(T, -2201514 // 234120865)
    a61 = convert(T, -7079 // 52152)
    a62 = convert(T, 767 // 2173)
    a63 = convert(T, 14027 // 52152)
    a64 = convert(T, 30 // 2173)
    #a65 = convert(T, 0 // 1)
    a71 = convert(T, 4817 // 51600)
    #a72 = convert(T, 0 // 1)
    a73 = convert(T, 388869 // 1216880)
    a74 = convert(T, 3276 // 23575)
    a75 = convert(T, -1142053 // 22015140)
    #a76 = convert(T, 0 // 1)
    abar21 = convert(T, 8 // 39)
    abar31 = convert(T, 1 // 13)
    abar32 = convert(T, 3 // 13)
    abar41 = convert(T, 7385 // 6912)
    abar42 = convert(T, -9425 // 2304)
    abar43 = convert(T, 13325 // 3456)
    abar51 = convert(T, 223324757 // 91364240)
    abar52 = convert(T, -174255393 // 18272848)
    abar53 = convert(T, 382840094 // 46824173)
    abar54 = convert(T, -39627252 // 234120865)
    abar61 = convert(T, 108475 // 36464)
    abar62 = convert(T, -9633 // 848)
    abar63 = convert(T, 7624604 // 806183)
    abar64 = convert(T, 8100 // 49979)
    abar65 = convert(T, -4568212 // 19446707)
    abar71 = convert(T, 4817 // 51600)
    #abar72 = convert(T, 0 // 1)
    abar73 = convert(T, 1685099 // 3650640)
    abar74 = convert(T, 19656 // 23575)
    abar75 = convert(T, -53676491 // 88060560)
    abar76 = convert(T, 53 // 240)
    b1 = convert(T, 4817 // 51600)
    #b2 = convert(T, 0 // 1)
    b3 = convert(T, 388869 // 1216880)
    b4 = convert(T, 3276 // 23575)
    b5 = convert(T, -1142053 // 22015140)
    #b6 = convert(T, 0 // 1)
    #b7 = convert(T, 0 // 1)
    bbar1 = convert(T, 4817 // 51600)
    #bbar2 = convert(T, 0 // 1)
    bbar3 = convert(T, 1685099 // 3650640)
    bbar4 = convert(T, 19656 // 23575)
    bbar5 = convert(T, -53676491 // 88060560)
    bbar6 = convert(T, 53 // 240)
    #bbar7 = convert(T, 0 // 1)
    Nystrom5ConstantCache(c1, c2, c3, c4, c5, c6, c7, a21, a31, a32, a41, a43, a51,
                          a52, a53, a54, a61, a62, a63, a64, a71, a73, a74, a75,
                          abar21, abar31, abar32, abar41, abar42, abar43, abar51,
                          abar52, abar53, abar54, abar61, abar62, abar63, abar64, abar65,
                          abar71, abar73, abar74, abar75, abar76, b1, b3, b4,
                          b5, bbar1, bbar3, bbar4, bbar5, bbar6)
end

struct IRKN3ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    bconst1::T
    bconst2::T
    c1::T2
    a21::T
    b1::T
    b2::T
    bbar1::T
    bbar2::T
end

function IRKN3ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    bconst1 = convert(T, 1.5)
    bconst2 = convert(T, -0.5)
    c1 = convert(T2, 0.5)
    a21 = convert(T, 0.125)
    b1 = convert(T, 0.6666666666666666)
    b2 = convert(T, 0.8333333333333334)
    bbar1 = convert(T, 0.3333333333333333)
    bbar2 = convert(T, 0.4166666666666667)
    IRKN3ConstantCache(bconst1, bconst2, c1, a21, b1, b2, bbar1, bbar2)
end

function IRKN3ConstantCache(T::Type, T2::Type)
    bconst1 = convert(T, 3 // 2)
    bconst2 = convert(T, -1 // 2)
    c1 = convert(T2, 1 // 2)
    a21 = convert(T, 1 // 8)
    b1 = convert(T, 2 // 3)
    b2 = convert(T, 5 // 6)
    bbar1 = convert(T, 1 // 3)
    bbar2 = convert(T, 5 // 12)
    IRKN3ConstantCache(bconst1, bconst2, c1, a21, b1, b2, bbar1, bbar2)
end

struct IRKN4ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    bconst1::T
    bconst2::T
    c1::T2
    c2::T2
    a21::T
    # a31::T
    a32::T
    b1::T
    b2::T
    b3::T
    bbar1::T
    bbar2::T
    bbar3::T
end

function IRKN4ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    bconst1 = convert(T, 1.5)
    bconst2 = convert(T, -0.5)
    c1 = convert(T2, 0.25)
    c2 = convert(T2, 0.75)
    a21 = convert(T, 0.03125)
    # a31     = convert(T,0)
    a32 = convert(T, 0.28125)
    b1 = convert(T, 1.0555555555555556)
    b2 = convert(T, -0.16666666666666666)
    b3 = convert(T, 0.6111111111111112)
    bbar1 = convert(T, -0.05555555555555555)
    bbar2 = convert(T, 0.2916666666666667)
    bbar3 = convert(T, 0.125)
    IRKN4ConstantCache(bconst1, bconst2, c1, c2, a21, a32, b1, b2, b3, bbar1, bbar2, bbar3)
end

function IRKN4ConstantCache(T::Type, T2::Type)
    bconst1 = convert(T, 3 // 2)
    bconst2 = convert(T, -1 // 2)
    c1 = convert(T2, 1 // 4)
    c2 = convert(T2, 3 // 4)
    a21 = convert(T, 1 // 32)
    # a31     = convert(T,0)
    a32 = convert(T, 9 // 32)
    b1 = convert(T, 19 // 18)
    b2 = convert(T, -1 // 6)
    b3 = convert(T, 11 // 18)
    bbar1 = convert(T, -1 // 18)
    bbar2 = convert(T, 7 // 24)
    bbar3 = convert(T, 1 // 8)
    IRKN4ConstantCache(bconst1, bconst2, c1, c2, a21, a32, b1, b2, b3, bbar1, bbar2, bbar3)
end

struct Nystrom5VelocityIndependentConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    c1::T2
    c2::T2
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    bbar1::T
    bbar2::T
    bbar3::T
    b1::T
    b2::T
    b3::T
    b4::T
end

function Nystrom5VelocityIndependentConstantCache(T::Type{<:CompiledFloats},
                                                  T2::Type{<:CompiledFloats})
    c1 = convert(T2, 0.2)
    c2 = convert(T2, 0.6666666666666666)
    # c3    = convert(T2,1)
    a21 = convert(T, 0.02)
    a31 = convert(T, -0.037037037037037035)
    a32 = convert(T, 0.25925925925925924)
    a41 = convert(T, 0.3)
    a42 = convert(T, -0.05714285714285714)
    a43 = convert(T, 0.2571428571428571)
    bbar1 = convert(T, 0.041666666666666664)
    bbar2 = convert(T, 0.2976190476190476)
    bbar3 = convert(T, 0.16071428571428573)
    b1 = bbar1
    b2 = convert(T, 0.37202380952380953)
    b3 = convert(T, 0.48214285714285715)
    b4 = convert(T, 0.10416666666666667)
    Nystrom5VelocityIndependentConstantCache(c1, c2, a21, a31, a32, a41, a42, a43, bbar1,
                                             bbar2, bbar3, b1, b2, b3, b4)
end

function Nystrom5VelocityIndependentConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 1 // 5)
    c2 = convert(T2, 2 // 3)
    # c3    = convert(T2,1)
    a21 = convert(T, 1 // 50)
    a31 = convert(T, -1 // 27)
    a32 = convert(T, 7 // 27)
    a41 = convert(T, 3 // 10)
    a42 = convert(T, -2 // 35)
    a43 = convert(T, 9 // 35)
    bbar1 = convert(T, 14 // 336)
    bbar2 = convert(T, 100 // 336)
    bbar3 = convert(T, 54 // 336)
    b1 = bbar1
    b2 = convert(T, 125 // 336)
    b3 = convert(T, 162 // 336)
    b4 = convert(T, 35 // 336)
    Nystrom5VelocityIndependentConstantCache(c1, c2, a21, a31, a32, a41, a42, a43, bbar1,
                                             bbar2, bbar3, b1, b2, b3, b4)
end

struct ERKN4ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    c1::T2
    c2::T2
    c3::T2
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    b1::T
    b2::T
    b3::T
    b4::T
    bp1::T # bp denotes bprime
    bp2::T
    bp3::T
    bp4::T
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    bptilde1::T
    bptilde2::T
    bptilde3::T
    bptilde4::T
end

function ERKN4ConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 1 // 4)
    c2 = convert(T2, 7 // 10)
    c3 = convert(T2, 1)
    a21 = convert(T, 1 // 32)
    a31 = convert(T, 19 // 600)
    a32 = convert(T, 16 // 75)
    a41 = convert(T, 32 // 315)
    a42 = convert(T, 58 // 315)
    a43 = convert(T, 3 // 14)
    btilde1 = convert(T, 1 // 21 - 14 // 375)
    btilde2 = convert(T, 28 // 81 - 136 // 375)
    btilde3 = convert(T, 50 // 567 - 2 // 25)
    btilde4 = convert(T, 1 // 54 - 1 // 50)
    bptilde1 = convert(T, 1 // 14 - 17 // 231)
    bptilde2 = convert(T, 32 // 81 - 116 // 297)
    bptilde3 = convert(T, 250 // 567 - 925 // 2079)
    bptilde4 = convert(T, 5 // 54 - 1 // 11)
    b1 = convert(T, 1 // 21)
    b2 = convert(T, 28 // 81)
    b3 = convert(T, 50 // 567)
    b4 = convert(T, 1 // 54)
    bp1 = convert(T, 1 // 14)
    bp2 = convert(T, 32 // 81)
    bp3 = convert(T, 250 // 567)
    bp4 = convert(T, 5 // 54)
    ERKN4ConstantCache(c1, c2, c3, a21, a31, a32, a41, a42, a43, b1, b2, b3, b4, bp1, bp2,
                       bp3, bp4, btilde1, btilde2, btilde3, btilde4, bptilde1, bptilde2,
                       bptilde3, bptilde4)
end

function ERKN4ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    ERKN4ConstantCache(convert(T2, 0.25),
                       convert(T2, 0.7),
                       convert(T2, 1.0),
                       convert(T, 0.03125),
                       convert(T, 0.03166666666666667),
                       convert(T, 0.21333333333333335),
                       convert(T, 0.10158730158730159),
                       convert(T, 0.18412698412698414),
                       convert(T, 0.21428571428571427),
                       convert(T, 0.047619047619047616),
                       convert(T, 0.345679012345679),
                       convert(T, 0.08818342151675485),
                       convert(T, 0.018518518518518517),
                       convert(T, 0.07142857142857142),
                       convert(T, 0.3950617283950617),
                       convert(T, 0.4409171075837742),
                       convert(T, 0.09259259259259259),
                       convert(T, 0.010285714285714285),
                       convert(T, -0.016987654320987654),
                       convert(T, 0.00818342151675485),
                       convert(T, -0.0014814814814814814),
                       convert(T, -0.0021645021645021645),
                       convert(T, 0.004489337822671156),
                       convert(T, -0.004008337341670675),
                       convert(T, 0.0016835016835016834))
end

struct ERKN5ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    c1::T2
    c2::T2
    c3::T2
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    b1::T
    b2::T
    b3::T
    b4::T
    bp1::T # bp denotes bprime
    bp2::T
    bp3::T
    bp4::T
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    # bptilde1::T
    # bptilde2::T
    # bptilde3::T
    # bptilde4::T
end

function ERKN5ConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 1 // 2)
    c2 = convert(T2, 19 // 70)
    c3 = convert(T2, 44 // 51)
    a21 = convert(T, 1 // 8)
    a31 = convert(T, 2907 // 343000)
    a32 = convert(T, 1216 // 42875)
    a41 = convert(T, 6624772 // Int64(128538819))
    a42 = convert(T, 6273905 // Int64(54121608))
    a43 = convert(T, Int64(210498365) // Int64(1028310552))
    b1 = convert(T, 479 // 5016)
    b2 = convert(T, 235 // 1776)
    b3 = convert(T, 145775 // 641744)
    b4 = convert(T, 309519 // 6873416)
    btilde1 = convert(T, 479 // 5016 - 184883 // 2021250)
    btilde2 = convert(T, 235 // 1776 - 411163 // 3399375)
    btilde3 = convert(T, 145775 // 641744 - 6 // 25)
    btilde4 = convert(T, 309519 // 6873416 - 593028 // Int64(12464375))
    bp1 = b1
    bp2 = convert(T, 235 // 888)
    bp3 = convert(T, 300125 // 962616)
    bp4 = convert(T, 2255067 // 6873416)
    # bptilde1 = convert(T,0)
    # bptilde2 = convert(T,0)
    # bptilde3 = convert(T,0)
    # bptilde4 = convert(T,0)
    ERKN5ConstantCache(c1, c2, c3, a21, a31, a32, a41, a42, a43, b1, b2, b3, b4, bp1, bp2,
                       bp3, bp4, btilde1, btilde2, btilde3, btilde4)
end

function ERKN5ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    ERKN5ConstantCache(convert(T2, 0.5),
                       convert(T2, 0.2714285714285714),
                       convert(T2, 0.8627450980392157),
                       convert(T, 0.125),
                       convert(T, 0.008475218658892128),
                       convert(T, 0.028361516034985424),
                       convert(T, 0.051539076300366506),
                       convert(T, 0.11592236875149756),
                       convert(T, 0.20470310704348388),
                       convert(T, 0.09549441786283891),
                       convert(T, 0.13231981981981983),
                       convert(T, 0.22715444164651324),
                       convert(T, 0.04503132067082801),
                       convert(T, 0.09549441786283891),
                       convert(T, 0.26463963963963966),
                       convert(T, 0.3117806061814888),
                       convert(T, 0.32808533631603265),
                       convert(T, 0.004024782736060931),
                       convert(T, 0.011367291781577495),
                       convert(T, -0.012845558353486749),
                       convert(T, -0.0025465161641516788))
end

struct ERKN7ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    a51::T
    a52::T
    a53::T
    a54::T
    a61::T
    a62::T
    a63::T
    a64::T
    a65::T
    a71::T
    a73::T
    a74::T
    a75::T
    a76::T
    b1::T
    b3::T
    b4::T
    b5::T
    b6::T
    bp1::T # bp denotes bprime
    bp3::T
    bp4::T
    bp5::T
    bp6::T
    bp7::T
    btilde1::T
    btilde3::T
    btilde4::T
    btilde5::T
    btilde6::T
    bptilde1::T
    bptilde3::T
    bptilde4::T
    bptilde5::T
    bptilde6::T
    bptilde7::T
end

function ERKN7ConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 108816483 // 943181462)
    c2 = convert(T2, 108816483 // 471590731)
    c3 = convert(T2, 151401202 // 200292705)
    c4 = convert(T2, 682035803 // 631524599)
    c5 = convert(T2, 493263404 // 781610081)
    c6 = convert(T2, 1)
    a21 = convert(T, 5107771 // 767472028)
    a31 = convert(T, 5107771 // 575604021)
    a32 = convert(T, 16661485 // 938806552)
    a41 = convert(T, 325996677 // 876867260)
    a42 = convert(T, -397622579 // 499461366)
    a43 = convert(T, 541212017 // 762248206)
    a51 = convert(T, 82243160 // 364375691)
    a52 = convert(T, -515873404 // 1213273815)
    a53 = convert(T, 820109726 // 1294837243)
    a54 = convert(T, 36245507 // 242779260)
    a61 = convert(T, 3579594 // 351273191)
    a62 = convert(T, 34292133 // 461028419)
    a63 = convert(T, 267156948 // 2671391749)
    a64 = convert(T, 22665163 // 1338599875)
    a65 = convert(T, -3836509 // 1614789462)
    a71 = convert(T, 53103334 // 780726093)
    a73 = convert(T, 352190060 // 1283966121)
    a74 = convert(T, 37088117 // 2206150964)
    a75 = convert(T, 7183323 // 1828127386)
    a76 = convert(T, 187705681 // 1370684829)
    b1 = convert(T, 53103334 // 780726093)
    b3 = convert(T, 352190060 // 1283966121)
    b4 = convert(T, 37088117 // 2206150964)
    b5 = convert(T, 7183323 // 1828127386)
    b6 = convert(T, 187705681 // 1370684829)
    bp1 = convert(T, 53103334 // 780726093)
    bp3 = convert(T, 244481296 // 685635505)
    bp4 = convert(T, 41493456 // 602487871)
    bp5 = convert(T, -45498718 // 926142189)
    bp6 = convert(T, 1625563237 // 4379140271)
    bp7 = convert(T, 191595797 // 1038702495)
    btilde1 = convert(T, 53103334 // 780726093 - 41808761 // 935030896)
    btilde3 = convert(T, 352190060 // 1283966121 - 46261019 // 135447428)
    btilde4 = convert(T, 37088117 // 2206150964 - 289298425 // 1527932372)
    btilde5 = convert(T, 7183323 // 1828127386 + 52260067 // 3104571287)
    btilde6 = convert(T, 187705681 // 1370684829 + 49872919 // 848719175)
    bptilde1 = convert(T, 53103334 // 780726093 - 41808761 // 935030896)
    bptilde3 = convert(T, 244481296 // 685635505 - 224724272 // 506147085)
    bptilde4 = convert(T, 41493456 // 602487871 - 2995752066 // 3862177123)
    bptilde5 = convert(T, -45498718 // 926142189 - 170795979 // 811534085)
    bptilde6 = convert(T, 1625563237 // 4379140271 + 177906423 // 1116903503)
    bptilde7 = convert(T, 191595797 // 1038702495 + 655510901 // 2077404990)
    ERKN7ConstantCache(c1, c2, c3, c4, c5, c6, a21, a31, a32, a41, a42, a43, a51, a52, a53,
                       a54, a61, a62, a63, a64, a65, a71, a73, a74, a75, a76, b1, b3, b4,
                       b5,
                       b6, bp1, bp3, bp4, bp5, bp6, bp7, btilde1, btilde3, btilde4, btilde5,
                       btilde6, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6, bptilde7)
end

function ERKN7ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    ERKN7ConstantCache(convert(T2, 108816483 // 943181462),
                       convert(T2, 0.23074347277618568),
                       convert(T2, 0.7558997318449516),
                       convert(T2, 1.0799829556599743),
                       convert(T2, 0.6310862871278652),
                       convert(T2, 1.0),
                       convert(T, 0.006655318778601792),
                       convert(T, 0.008873758371469056),
                       convert(T, 0.01774751674293811),
                       convert(T, 0.37177426033673555),
                       convert(T, -0.796102774043188),
                       convert(T, 0.7100207160080872),
                       convert(T, 0.2257097880879216),
                       convert(T, -0.4251912450611983),
                       convert(T, 0.6333689662029593),
                       convert(T, 0.14929408302834435),
                       convert(T, 0.010190342137439119),
                       convert(T, 0.07438182026691938),
                       convert(T, 0.10000665312379087),
                       convert(T, 0.016931992467129134),
                       convert(T, -0.002375857094861324),
                       convert(T, 0.06801788037587723),
                       convert(T, 0.2742985614960786),
                       convert(T, 0.01681123259704543),
                       convert(T, 0.003929333948504177),
                       convert(T, 0.13694299158249457),
                       convert(T, 0.06801788037587723),
                       convert(T, 0.2742985614960786),
                       convert(T, 0.01681123259704543),
                       convert(T, 0.003929333948504177),
                       convert(T, 0.13694299158249457),
                       convert(T, 0.06801788037587723),
                       convert(T, 0.35657618985177847),
                       convert(T, 0.06887019307314819),
                       convert(T, -0.049127141102520276),
                       convert(T, 0.371206021365649),
                       convert(T, 0.18445685643606738),
                       convert(T, 0.023304105484742516),
                       convert(T, -0.06724368617214582),
                       convert(T, -0.1725285773981577),
                       convert(T, 0.020762597600343137),
                       convert(T, 0.1957055604852179),
                       convert(T, 0.023304105484742516),
                       convert(T, -0.08741386493634701),
                       convert(T, -0.7067938872093759),
                       convert(T, -0.25958777628335805),
                       convert(T, 0.5304914229443385),
                       convert(T, 0.5))
end

struct DPRKN4ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    c1::T2
    c2::T2
    c3::T2
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    b1::T
    b2::T
    b3::T
    bp1::T # bp denotes bprime
    bp2::T
    bp3::T
    bp4::T
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    bptilde1::T
    bptilde2::T
    bptilde3::T
    bptilde4::T
end

function DPRKN4ConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 1 // 4)
    c2 = convert(T2, 7 // 10)
    c3 = convert(T2, 1)
    a21 = convert(T, 1 // 32)
    a31 = convert(T, 7 // 1000)
    a32 = convert(T, 119 // 500)
    a41 = convert(T, 1 // 14)
    a42 = convert(T, 8 // 27)
    a43 = convert(T, 25 // 189)
    b1 = convert(T, 1 // 14)
    b2 = convert(T, 8 // 27)
    b3 = convert(T, 25 // 189)
    # b4 = convert(T, 0)
    bp1 = convert(T, 1 // 14)
    bp2 = convert(T, 32 // 81)
    bp3 = convert(T, 250 // 567)
    bp4 = convert(T, 5 // 54)
    btilde1 = convert(T, 1 // 14 + 7 // 150)
    btilde2 = convert(T, 8 // 27 - 67 // 150)
    btilde3 = convert(T, 25 // 189 - 3 // 20)
    btilde4 = convert(T, 1 // 20)
    bptilde1 = convert(T, 1 // 14 - 13 // 21)
    bptilde2 = convert(T, 32 // 81 + 20 // 27)
    bptilde3 = convert(T, 250 // 567 - 275 // 189)
    bptilde4 = convert(T, 5 // 54 + 1 // 3)
    DPRKN4ConstantCache(c1, c2, c3, a21, a31, a32, a41, a42, a43, b1, b2, b3,
                        bp1, bp2, bp3, bp4, btilde1, btilde2, btilde3, btilde4,
                        bptilde1, bptilde2, bptilde3, bptilde4)
end

function DPRKN4ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    c1 = convert(T2, 0.25)
    c2 = convert(T2, 0.7)
    c3 = convert(T2, 1.0)
    a21 = convert(T, 0.03125)
    a31 = convert(T, 0.007)
    a32 = convert(T, 0.238)
    a41 = convert(T, 0.07142857142857142)
    a42 = convert(T, 0.2962962962962963)
    a43 = convert(T, 0.13227513227513227)
    b1 = convert(T, 0.07142857142857142)
    b2 = convert(T, 0.2962962962962963)
    b3 = convert(T, 0.13227513227513227)
    bp1 = convert(T, 0.07142857142857142)
    bp2 = convert(T, 0.3950617283950617)
    bp3 = convert(T, 0.4409171075837742)
    bp4 = convert(T, 0.09259259259259259)
    btilde1 = convert(T, 0.11809523809523809)
    btilde2 = convert(T, -0.15037037037037038)
    btilde3 = convert(T, -0.017724867724867727)
    btilde4 = convert(T, 0.05)
    bptilde1 = convert(T, -0.5476190476190477)
    bptilde2 = convert(T, 1.1358024691358024)
    bptilde3 = convert(T, -1.0141093474426808)
    bptilde4 = convert(T, 0.42592592592592593)
    DPRKN4ConstantCache(c1, c2, c3, a21, a31, a32, a41, a42, a43, b1, b2, b3,
                        bp1, bp2, bp3, bp4, btilde1, btilde2, btilde3, btilde4,
                        bptilde1, bptilde2, bptilde3, bptilde4)
end
struct DPRKN5ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    a21::T
    a31::T
    a32::T
    a41::T
    # a42::T
    a43::T
    a51::T
    # a52::T
    a53::T
    a54::T
    a61::T
    # a62::T
    a63::T
    a64::T
    a65::T
    b1::T
    # b2::T
    b3::T
    b4::T
    b5::T
    # b6::T
    bp1::T # bp denotes bprime
    # bp2::T
    bp3::T
    bp4::T
    bp5::T
    bp6::T
    btilde1::T
    # btilde2::T
    btilde3::T
    btilde4::T
    btilde5::T
    # btilde6::T
    bptilde1::T
    # bptilde2::T
    bptilde3::T
    bptilde4::T
    bptilde5::T
    bptilde6::T
end

function DPRKN5ConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 1 // 8)
    c2 = convert(T2, 1 // 4)
    c3 = convert(T2, 1 // 2)
    c4 = convert(T2, 3 // 4)
    c5 = convert(T2, 1)
    a21 = convert(T, 1 // 128)
    a31 = convert(T, 1 // 96)
    a32 = convert(T, 1 // 48)
    a41 = convert(T, 1 // 24)
    # a42 = convert(T, 0)
    a43 = convert(T, 1 // 12)
    a51 = convert(T, 9 // 128)
    # a52 = convert(T, 0)
    a53 = convert(T, 9 // 64)
    a54 = convert(T, 9 // 128)
    a61 = convert(T, 7 // 90)
    # a62 = convert(T, 0)
    a63 = convert(T, 4 // 15)
    a64 = convert(T, 1 // 15)
    a65 = convert(T, 4 // 45)
    b1 = convert(T, 7 // 90)
    # b2 = convert(T,0)
    b3 = convert(T, 4 // 15)
    b4 = convert(T, 1 // 15)
    b5 = convert(T, 4 // 45)
    # b6 = convert(T, 0)
    bp1 = convert(T, 7 // 90)
    # bp2 = convert(T,0)
    bp3 = convert(T, 16 // 45)
    bp4 = convert(T, 2 // 15)
    bp5 = convert(T, 16 // 45)
    bp6 = convert(T, 7 // 90)
    btilde1 = convert(T, 7 // 90 - 1 // 6)
    # btilde2 = convert(T,0)
    btilde3 = convert(T, 4 // 15)
    btilde4 = convert(T, 1 // 15 - 1 // 3)
    btilde5 = convert(T, 4 // 45)
    #btilde6 = convert(T, 0)
    bptilde1 = convert(T, 7 // 90)
    # bptilde2 = convert(T,0)
    bptilde3 = convert(T, 16 // 45 - 2 // 3)
    bptilde4 = convert(T, 2 // 15 + 1 // 3)
    bptilde5 = convert(T, 16 // 45 - 2 // 3)
    bptilde6 = convert(T, 7 // 90)
    DPRKN5ConstantCache(c1, c2, c3, c4, c5, a21, a31, a32, a41, a43, a51,
                        a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1,
                        bp3, bp4, bp5, bp6, btilde1, btilde3, btilde4, btilde5,
                        bptilde1, bptilde3, bptilde4, bptilde5, bptilde6)
end

function DPRKN5ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    c1 = convert(T2, 0.125)
    c2 = convert(T2, 0.25)
    c3 = convert(T2, 0.5)
    c4 = convert(T2, 0.75)
    c5 = convert(T2, 1.0)
    a21 = convert(T, 1 // 128)
    a31 = convert(T, 1 // 96)
    a32 = convert(T, 1 // 48)
    a41 = convert(T, 1 // 24)
    a43 = convert(T, 1 // 12)
    a51 = convert(T, 7 // 90)
    a53 = convert(T, 4 // 15)
    a54 = convert(T, 1 // 15)
    a61 = convert(T, 0.07777777777777778)
    a63 = convert(T, 0.26666666666666666)
    a64 = convert(T, 0.06666666666666667)
    a65 = convert(T, 0.08888888888888889)
    b1 = convert(T, 0.07777777777777778)
    b3 = convert(T, 0.26666666666666666)
    b4 = convert(T, 0.06666666666666667)
    b5 = convert(T, 0.08888888888888889)
    bp1 = convert(T, 0.07777777777777778)
    bp3 = convert(T, 0.35555555555555557)
    bp4 = convert(T, 0.13333333333333333)
    bp5 = convert(T, 0.35555555555555557)
    bp6 = convert(T, 0.07777777777777778)
    btilde1 = convert(T, -0.08888888888888888)
    btilde3 = convert(T, 0.26666666666666666)
    btilde4 = convert(T, -0.26666666666666666)
    btilde5 = convert(T, 0.08888888888888889)
    bptilde1 = convert(T, 0.07777777777777778)
    bptilde3 = convert(T, -0.31111111111111106)
    bptilde4 = convert(T, 0.4666666666666667)
    bptilde5 = convert(T, -0.31111111111111106)
    bptilde6 = convert(T, 0.07777777777777778)
    DPRKN5ConstantCache(c1, c2, c3, c4, c5, a21, a31, a32, a41, a43, a51,
                        a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1,
                        bp3, bp4, bp5, bp6, btilde1, btilde3, btilde4, btilde5,
                        bptilde1, bptilde3, bptilde4, bptilde5, bptilde6)
end

struct DPRKN6ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    a51::T
    a52::T
    a53::T
    a54::T
    a61::T
    # a62::T
    a63::T
    a64::T
    a65::T
    b1::T
    # b2::T
    b3::T
    b4::T
    b5::T
    # b6::T
    bp1::T # bp denotes bprime
    # bp2::T
    bp3::T
    bp4::T
    bp5::T
    bp6::T
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    btilde5::T
    # btilde6::T
    bptilde1::T
    # bptilde2::T
    bptilde3::T
    bptilde4::T
    bptilde5::T
    bptilde6::T
    r14::T
    r13::T
    r12::T
    r11::T
    r10::T
    r34::T
    r33::T
    r32::T
    r31::T
    r44::T
    r43::T
    r42::T
    r41::T
    r54::T
    r53::T
    r52::T
    r51::T
    r64::T
    r63::T
    r62::T
    r61::T
    rp14::T
    rp13::T
    rp12::T
    rp11::T
    rp10::T
    rp34::T
    rp33::T
    rp32::T
    rp31::T
    rp44::T
    rp43::T
    rp42::T
    rp41::T
    rp54::T
    rp53::T
    rp52::T
    rp51::T
    rp64::T
    rp63::T
    rp62::T
    rp61::T
end

function DPRKN6ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    c1 = convert(T2, 0.12929590313670442)
    c2 = convert(T2, 0.25859180627340883)
    c3 = convert(T2, 0.67029708261548)
    c4 = convert(T2, 0.9)
    c5 = convert(T2, 1.0)
    a21 = convert(T, 0.008358715283968025)
    a31 = convert(T, 0.011144953711957367)
    a32 = convert(T, 0.022289907423914734)
    a41 = convert(T, 0.1454747428010918)
    a42 = convert(T, -0.22986064052264749)
    a43 = convert(T, 0.3090349872029675)
    a51 = convert(T, -0.20766826295078997)
    a52 = convert(T, 0.6863667842925143)
    a53 = convert(T, -0.19954927787234925)
    a54 = convert(T, 0.12585075653062489)
    a61 = convert(T, 0.07811016144349478)
    a63 = convert(T, 0.2882917411897668)
    a64 = convert(T, 0.12242553717457041)
    a65 = convert(T, 0.011172560192168035)
    b1 = convert(T, 0.07811016144349478)
    b3 = convert(T, 0.2882917411897668)
    b4 = convert(T, 0.12242553717457041)
    b5 = convert(T, 0.011172560192168035)
    bp1 = convert(T, 0.07811016144349478)
    bp3 = convert(T, 0.3888434787059826)
    bp4 = convert(T, 0.3713207579288423)
    bp5 = convert(T, 0.11172560192168035)
    bp6 = convert(T, 0.05)
    btilde1 = convert(T, -0.9807490989269235)
    btilde2 = convert(T, 2.406751371924452)
    btilde3 = convert(T, -1.559600370364267)
    btilde4 = convert(T, 0.12242553717457041)
    btilde5 = convert(T, 0.011172560192168035)
    bptilde1 = convert(T, 0.023504273504273504)
    bptilde3 = convert(T, -0.07242330719764424)
    bptilde4 = convert(T, 0.17543989844952962)
    bptilde5 = convert(T, -0.2765208647561589)
    bptilde6 = convert(T, 0.15)
    r14 = convert(T, 0.21367521367521367)
    r13 = convert(T, -0.9066951566951567)
    r12 = convert(T, 1.5161443494776827)
    r11 = convert(T, -1.245014245014245)
    r10 = convert(T, 0.5)
    r34 = convert(T, -0.6583937017967658)
    r33 = convert(T, 2.5384011164109506)
    r32 = convert(T, -3.577652872294921)
    r31 = convert(T, 1.9859371988705032)
    r44 = convert(T, 1.5949081677229964)
    r43 = convert(T, -5.164133553908094)
    r42 = convert(T, 5.547586751052329)
    r41 = convert(T, -1.8559358276926614)
    r54 = convert(T, -2.513826043237808)
    r53 = convert(T, 7.273336685101391)
    r52 = convert(T, -6.926987319144182)
    r51 = convert(T, 2.178649237472767)
    r64 = convert(T, 1.3636363636363635)
    r63 = convert(T, -3.7409090909090907)
    r62 = convert(T, 3.440909090909091)
    r61 = convert(T, -1.0636363636363637)
    rp14 = convert(T, 1.2820512820512822)
    rp13 = convert(T, -4.533475783475783)
    rp12 = convert(T, 6.064577397910731)
    rp11 = convert(T, -3.735042735042735)
    rp10 = convert(T, 1)
    rp34 = convert(T, -3.950362210780595)
    rp33 = convert(T, 12.692005582054751)
    rp32 = convert(T, -14.310611489179683)
    rp31 = convert(T, 5.95781159661151)
    rp44 = convert(T, 9.56944900633798)
    rp43 = convert(T, -25.820667769540467)
    rp42 = convert(T, 22.190347004209315)
    rp41 = convert(T, -5.567807483077984)
    rp54 = convert(T, -15.082956259426847)
    rp53 = convert(T, 36.366683425506956)
    rp52 = convert(T, -27.707949276576727)
    rp51 = convert(T, 6.5359477124183005)
    rp64 = convert(T, 8.181818181818182)
    rp63 = convert(T, -18.704545454545453)
    rp62 = convert(T, 13.763636363636364)
    rp61 = convert(T, -3.190909090909091)
    DPRKN6ConstantCache(c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51,
                        a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1,
                        bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4,
                        btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6,
                        r14, r13, r12, r11, r10, r34, r33, r32, r31, r44, r43, r42, r41,
                        r54,
                        r53, r52, r51, r64, r63, r62, r61, rp14, rp13, rp12, rp11, rp10,
                        rp34,
                        rp33, rp32, rp31, rp44, rp43, rp42, rp41, rp54, rp53, rp52, rp51,
                        rp64, rp63, rp62, rp61)
end

function DPRKN6ConstantCache(T::Type, T2::Type)
    R = sqrt(big(8581))
    c1 = convert(T2, (209 - R) / 900)
    c2 = convert(T2, (209 - R) / 450)
    c3 = convert(T2, (209 + R) / 450)
    c4 = convert(T2, 9 // 10)
    c5 = convert(T2, 1)
    a21 = convert(T, (26131 - 209R) / 81_0000)
    a31 = convert(T, (26131 - 209R) / 60_7500)
    a32 = convert(T, (26131 - 209R) / 30_3750)
    a41 = convert(T, (980403512254 + 7781688431R) / 116944_6992_1875)
    a42 = convert(T, -(126288_4486208 + 153854_81287R) / 116944_6992_1875)
    a43 = convert(T, (7166_233_891_441 + 786_945_632_99R) / 46_777_879_687_500)
    a51 = convert(T, -9(329260 + 3181R) / 2704_0000)
    a52 = convert(T, 27(35129 + 3331R) / 1352_0000)
    a53 = convert(T, -27(554358343 + 31040327R) / 46406048_0000)
    a54 = convert(T, 153(8555_257 - 67973R) / 274592_0000)
    a61 = convert(T, 329 // 4212)
    # a62      = convert(T,0)
    a63 = convert(T, (8411_9543 + 366_727R) / 4096_22616)
    a64 = convert(T, (8411_9543 - 366_727R) / 4096_22616)
    a65 = convert(T, 200 // 17901)
    b1 = convert(T, 329 // 4212)
    # b2       = convert(T,0)
    b3 = a63
    b4 = a64
    b5 = convert(T, 200 // 17901)
    # b6       = convert(T,0)
    bp1 = b1
    # bp2      = b2
    bp3 = convert(T, (389225579 + 96856R) / 10_2405_6540)
    bp4 = convert(T, (389225579 - 96856R) / 10_2405_6540)
    bp5 = convert(T, 2000 // 17901)
    bp6 = convert(T, 1 // 20)
    btilde1 = convert(T, 329 // 4212 - (2701 + 23R) / 4563)
    btilde2 = convert(T, (9829 + 131R) / 9126)
    btilde3 = convert(T, (8411_9543 + 366_727R) / 4096_22616 - 5(1798 + 17R) / 9126)
    btilde4 = b4
    btilde5 = b5
    # btilde6  = convert(T,0)
    bptilde1 = convert(T, 329 // 4212 - 115 // 2106)
    # btildep2 = convert(T,0)
    bptilde3 = convert(T,
                       (389225579 + 96856R) / 10_2405_6540 -
                       (8411_9543 + 366_727R) / 2560_14135)
    bptilde4 = convert(T,
                       (389225579 - 96856R) / 10_2405_6540 -
                       (8411_9543 - 366_727R) / 2560_14135)
    bptilde5 = convert(T, 2000 // 17901 - 6950 // 17901)
    bptilde6 = convert(T, 1 // 20 + 1 // 10)
    r14 = convert(T, 900 // 4212)
    r13 = convert(T, -3819 // 4212)
    r12 = convert(T, 6386 // 4212)
    r11 = convert(T, -5244 // 4212)
    r10 = convert(T, 2106 // 4212)
    r34 = convert(T, 1800 * (5860823 - 152228R) / 22529243880)
    r33 = convert(T, -6 * (4929647204 - 156109769R) / 22529243880)
    r32 = convert(T, (22190560391 - 1109665151R) / 22529243880)
    r31 = convert(T, 18 * (81356461 + 25954829R) / 22529243880)
    r44 = convert(T, 1800 * (5860823 + 152228R) / 22529243880)
    r43 = convert(T, -6 * (4929647204 + 156109769R) / 22529243880)
    r42 = convert(T, (22190560391 + 1109665151R) / 22529243880)
    r41 = convert(T, 18 * (81356461 - 25954829R) / 22529243880)
    r54 = convert(T, -200 * 225 // 17901)
    r53 = convert(T, 200 * 651 // 17901)
    r52 = convert(T, -200 * 620 // 17901)
    r51 = convert(T, 200 * 195 // 17901)
    r64 = convert(T, 15 // 11)
    r63 = convert(T, -823 // 220)
    r62 = convert(T, 757 // 220)
    r61 = convert(T, -117 // 110)
    rp14 = convert(T, 5400 // 4212)
    rp13 = convert(T, -19095 // 4212)
    rp12 = convert(T, 25544 // 4212)
    rp11 = convert(T, -15732 // 4212)
    rp10 = convert(T, 1)
    rp34 = convert(T, 5400 * (5860823 - 152228R) / 11264621940)
    rp33 = convert(T, -15 * (4929647204 - 156109769R) / 11264621940)
    rp32 = convert(T, 2 * (22190560391 - 1109665151R) / 11264621940)
    rp31 = convert(T, 27 * (81356461 + 25954829R) / 11264621940)
    rp44 = convert(T, 5400 * (5860823 + 152228R) / 11264621940)
    rp43 = convert(T, -15 * (4929647204 + 156109769R) / 11264621940)
    rp42 = convert(T, 2 * (22190560391 + 1109665151R) / 11264621940)
    rp41 = convert(T, 27 * (81356461 - 25954829R) / 11264621940)
    rp54 = convert(T, -1000 * 270 // 17901)
    rp53 = convert(T, 1000 * 651 // 17901)
    rp52 = convert(T, -1000 * 496 // 17901)
    rp51 = convert(T, 1000 * 117 // 17901)
    rp64 = convert(T, 1800 // 220)
    rp63 = convert(T, -4115 // 220)
    rp62 = convert(T, 3028 // 220)
    rp61 = convert(T, -702 // 220)
    DPRKN6ConstantCache(c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51,
                        a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1,
                        bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4,
                        btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6,
                        r14, r13, r12, r11, r10, r34, r33, r32, r31, r44, r43, r42, r41,
                        r54,
                        r53, r52, r51, r64, r63, r62, r61, rp14, rp13, rp12, rp11, rp10,
                        rp34,
                        rp33, rp32, rp31, rp44, rp43, rp42, rp41, rp54, rp53, rp52, rp51,
                        rp64, rp63, rp62, rp61)
end

struct DPRKN6FMConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    a51::T
    a52::T
    a53::T
    a54::T
    a61::T
    a62::T
    a63::T
    a64::T
    a65::T
    b1::T
    b2::T
    b3::T
    b4::T
    b5::T
    # b6::T
    bp1::T # bp denotes bprime
    bp2::T
    bp3::T
    bp4::T
    bp5::T
    bp6::T
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    btilde5::T
    # btilde6::T
    bptilde1::T
    bptilde2::T
    bptilde3::T
    bptilde4::T
    bptilde5::T
    # bptilde6::T
end

function DPRKN6FMConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 1 // 10)
    c2 = convert(T2, 3 // 10)
    c3 = convert(T2, 7 // 10)
    c4 = convert(T2, 17 // 25)
    c5 = convert(T2, 1)
    a21 = convert(T, 1 // 200)
    a31 = convert(T, -1 // 2200)
    a32 = convert(T, 1 // 22)
    a41 = convert(T, 637 // 6600)
    a42 = convert(T, -7 // 110)
    a43 = convert(T, 7 // 33)
    a51 = convert(T, 225437 // 1968750)
    a52 = convert(T, -30073 // 281250)
    a53 = convert(T, 65569 // 281250)
    a54 = convert(T, -9367 // 984375)
    a61 = convert(T, 151 // 2142)
    a62 = convert(T, 5 // 116)
    a63 = convert(T, 385 // 1368)
    a64 = convert(T, 55 // 168)
    a65 = convert(T, -6250 // 28101)
    b1 = convert(T, 151 // 2142)
    b2 = convert(T, 5 // 116)
    b3 = convert(T, 385 // 1368)
    b4 = convert(T, 55 // 168)
    b5 = convert(T, -6250 // 28101)
    # b6 = convert(T, 0)
    bp1 = convert(T, 151 // 2142)
    bp2 = convert(T, 25 // 522)
    bp3 = convert(T, 275 // 684)
    bp4 = convert(T, 275 // 252)
    bp5 = convert(T, -78125 // 112404)
    bp6 = convert(T, 1 // 12)
    btilde1 = convert(T, 151 // 2142 - 1349 // 157500)
    btilde2 = convert(T, 5 // 116 - 7873 // 50000)
    btilde3 = convert(T, 385 // 1368 - 192199 // 900000)
    btilde4 = convert(T, 55 // 168 - 521683 // 2100000)
    btilde5 = convert(T, -6250 // 28101 + 16 // 125)
    # btilde6 = convert(T, 0)
    bptilde1 = convert(T, 151 // 2142 - 1349 // 157500)
    bptilde2 = convert(T, 25 // 522 - 7873 // 45000)
    bptilde3 = convert(T, 275 // 684 - 27457 // 90000)
    bptilde4 = convert(T, 275 // 252 - 521683 // 630000)
    bptilde5 = convert(T, -78125 // 112404 + 2 // 5)
    # bptilde6 = convert(T, 0)
    DPRKN6FMConstantCache(c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52,
                          a53, a54, a61, a62, a63, a64, a65, b1, b2, b3, b4, b5, bp1, bp2,
                          bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4, btilde5,
                          bptilde1, bptilde2, bptilde3, bptilde4, bptilde5)
end

function DPRKN6FMConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    c1 = convert(T2, 0.1)
    c2 = convert(T2, 0.3)
    c3 = convert(T2, 0.7)
    c4 = convert(T2, 0.68)
    c5 = convert(T2, 1.0)
    a21 = convert(T, 0.005)
    a31 = convert(T, -0.00045454545454545455)
    a32 = convert(T, 0.045454545454545456)
    a41 = convert(T, 0.09651515151515151)
    a42 = convert(T, -0.06363636363636363)
    a43 = convert(T, 0.21212121212121213)
    a51 = convert(T, 0.11450768253968253)
    a52 = convert(T, -0.10692622222222223)
    a53 = convert(T, 0.23313422222222221)
    a54 = convert(T, -0.00951568253968254)
    a61 = convert(T, 0.07049486461251167)
    a62 = convert(T, 0.04310344827586207)
    a63 = convert(T, 0.2814327485380117)
    a64 = convert(T, 0.3273809523809524)
    a65 = convert(T, -0.22241201380733783)
    b1 = convert(T, 0.07049486461251167)
    b2 = convert(T, 0.04310344827586207)
    b3 = convert(T, 0.2814327485380117)
    b4 = convert(T, 0.3273809523809524)
    b5 = convert(T, -0.22241201380733783)
    bp1 = convert(T, 0.07049486461251167)
    bp2 = convert(T, 0.04789272030651341)
    bp3 = convert(T, 0.402046783625731)
    bp4 = convert(T, 1.0912698412698412)
    bp5 = convert(T, -0.6950375431479306)
    bp6 = convert(T, 0.08333333333333333)
    btilde1 = convert(T, 0.061929785247432305)
    btilde2 = convert(T, -0.11435655172413792)
    btilde3 = convert(T, 0.06787830409356727)
    btilde4 = convert(T, 0.07896047619047619)
    btilde5 = convert(T, -0.09441201380733782)
    bptilde1 = convert(T, 0.061929785247432305)
    bptilde2 = convert(T, -0.12706283524904216)
    bptilde3 = convert(T, 0.0969690058479532)
    bptilde4 = convert(T, 0.26320158730158716)
    bptilde5 = convert(T, -0.2950375431479306)
    DPRKN6FMConstantCache(c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52,
                          a53, a54, a61, a62, a63, a64, a65, b1, b2, b3, b4, b5, bp1, bp2,
                          bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4, btilde5,
                          bptilde1, bptilde2, bptilde3, bptilde4, bptilde5)
end

struct DPRKN8ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2
    c7::T2
    c8::T2
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    a51::T
    a52::T
    a53::T
    a54::T
    a61::T
    a62::T
    a63::T
    a64::T
    a65::T
    a71::T
    a72::T
    a73::T
    a74::T
    a75::T
    a76::T
    a81::T
    a82::T
    a83::T
    a84::T
    a85::T
    a86::T
    a87::T
    a91::T
    # a92::T
    a93::T
    a94::T
    a95::T
    a96::T
    a97::T
    # a98::T
    b1::T
    # b2::T
    b3::T
    b4::T
    b5::T
    b6::T
    b7::T
    # b8::T
    # b9::T
    bp1::T
    # bp2::T
    bp3::T
    bp4::T
    bp5::T
    bp6::T
    bp7::T
    bp8::T
    # bp9::T
    btilde1::T
    # btilde2::T
    btilde3::T
    btilde4::T
    btilde5::T
    btilde6::T
    btilde7::T
    # btilde8::T
    # btilde9::T
    bptilde1::T
    # bptilde2::T
    bptilde3::T
    bptilde4::T
    bptilde5::T
    bptilde6::T
    bptilde7::T
    bptilde8::T
    bptilde9::T
end

function DPRKN8ConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 1 // 20)
    c2 = convert(T2, 1 // 10)
    c3 = convert(T2, 3 // 10)
    c4 = convert(T2, 1 // 2)
    c5 = convert(T2, 7 // 10)
    c6 = convert(T2, 9 // 10)
    c7 = convert(T2, 1)
    c8 = convert(T2, 1)
    a21 = convert(T, 1 // 800)
    a31 = convert(T, 1 // 600)
    a32 = convert(T, 1 // 300)
    a41 = convert(T, 9 // 200)
    a42 = convert(T, -9 // 100)
    a43 = convert(T, 9 // 100)
    a51 = convert(T, -66701 // 197352)
    a52 = convert(T, 28325 // 32892)
    a53 = convert(T, -2665 // 5482)
    a54 = convert(T, 2170 // 24669)
    a61 = convert(T, 2270_15747 // 30425_1000)
    a62 = convert(T, -5489_7451 // 30425_100)
    a63 = convert(T, 12942_349 // 10141_700)
    a64 = convert(T, -9499 // 304_251)
    a65 = convert(T, 539 // 9250)
    a71 = convert(T, -11318_91597 // 9017_89000)
    a72 = convert(T, 4196_4921 // 1288_2700)
    a73 = convert(T, -6663_147 // 3220_675)
    a74 = convert(T, 270_954 // 644_135)
    a75 = convert(T, -108 // 5875)
    a76 = convert(T, 114 // 1645)
    a81 = convert(T, 138_369_59 // 3667458)
    a82 = convert(T, -177_314_50 // 1833729)
    a83 = convert(T, 106_3919_505 // 15647_8208)
    a84 = convert(T, -332_138_45 // 3911_9552)
    a85 = convert(T, 133_35 // 285_44)
    a86 = convert(T, -705 // 14272)
    a87 = convert(T, 1645 // 57088)
    a91 = convert(T, 223 // 7938)
    # a92 = convert(T,0)
    a93 = convert(T, 1175 // 8064)
    a94 = convert(T, 925 // 6048)
    a95 = convert(T, 41 // 448)
    a96 = convert(T, 925 // 14112)
    a97 = convert(T, 1175 // 72576)
    # a98 = convert(T,0)
    b1 = convert(T, 223 // 7938)
    # b2 = convert(T,0)
    b3 = convert(T, 1175 // 8064)
    b4 = convert(T, 925 // 6048)
    b5 = convert(T, 41 // 448)
    b6 = convert(T, 925 // 14112)
    b7 = convert(T, 1175 // 72576)
    # b8 = convert(T,0)
    # b9 = convert(T,0)
    bp1 = convert(T, 223 // 7938)
    # bp2 = convert(T,0)
    bp3 = convert(T, 5875 // 36288)
    bp4 = convert(T, 4625 // 21168)
    bp5 = convert(T, 41 // 224)
    bp6 = convert(T, 4625 // 21168)
    bp7 = convert(T, 5875 // 36288)
    bp8 = convert(T, 223 // 7938)
    # bp9 = convert(T,0)
    btilde1 = convert(T, 223 // 7938 - 7987_313 // 10994_1300)
    # btilde2 = convert(T,0)
    btilde3 = convert(T, 1175 // 8064 - 1610_737 // 4467_4560)
    btilde4 = convert(T, 925 // 6048 - 10023_263 // 3350_5920)
    btilde5 = convert(T, 41 // 448 + 497_221 // 1240_9600)
    btilde6 = convert(T, 925 // 14112 - 1002_3263 // 7818_0480)
    btilde7 = convert(T, 1175 // 72576 - 1610_737 // 40207_1040)
    # btilde8 = convert(T,0)
    # btilde9 = convert(T,0)
    bptilde1 = convert(T, 223 // 7938 - 7987_313 // 10994_1300)
    # bptilde2 = convert(T,0)
    bptilde3 = convert(T, 5875 // 36288 - 1610_737 // 4020_7104)
    bptilde4 = convert(T, 4625 // 21168 - 1002_3263 // 2345_4144)
    bptilde5 = convert(T, 41 // 224 + 497_221 // 620_4800)
    bptilde6 = convert(T, 4625 // 21168 - 1002_3263 // 2345_4144)
    bptilde7 = convert(T, 5875 // 36288 - 1610_737 // 40207_104)
    bptilde8 = convert(T, 223 // 7938 + 4251_941 // 5497_0650)
    bptilde9 = convert(T, -3 // 20)
    DPRKN8ConstantCache(c1, c2, c3, c4, c5, c6, c7, c8, a21, a31, a32, a41, a42, a43, a51,
                        a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75,
                        a76, a81, a82, a83, a84, a85, a86, a87, a91, a93, a94, a95, a96,
                        a97, b1, b3, b4, b5, b6, b7, bp1, bp3, bp4, bp5, bp6, bp7, bp8,
                        btilde1, btilde3, btilde4, btilde5, btilde6, btilde7, bptilde1,
                        bptilde3, bptilde4, bptilde5, bptilde6, bptilde7, bptilde8,
                        bptilde9)
end

function DPRKN8ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    DPRKN8ConstantCache(convert(T2, 0.05),
                        convert(T2, 0.1),
                        convert(T2, 0.3),
                        convert(T2, 0.5),
                        convert(T2, 0.7),
                        convert(T2, 0.9),
                        convert(T2, 1.0),
                        convert(T2, 1.0),
                        convert(T, 0.00125),
                        convert(T, 0.0016666666666666668),
                        convert(T, 0.0033333333333333335),
                        convert(T, 0.045),
                        convert(T, -0.09),
                        convert(T, 0.09),
                        convert(T, -0.3379798532571243),
                        convert(T, 0.8611516478170984),
                        convert(T, -0.48613644655235316),
                        convert(T, 0.0879646519923791),
                        convert(T, 0.7461462641043086),
                        convert(T, -1.804347430246737),
                        convert(T, 1.2761518285888953),
                        convert(T, -0.031220932716737166),
                        convert(T, 0.05827027027027027),
                        convert(T, -1.2551623461807584),
                        convert(T, 3.257463187064823),
                        convert(T, -2.068866619575089),
                        convert(T, 0.4206478455603251),
                        convert(T, -0.018382978723404254),
                        convert(T, 0.06930091185410335),
                        convert(T, 3.772901830095941),
                        convert(T, -9.669613121677195),
                        convert(T, 6.7991544547851674),
                        convert(T, -0.8490343907823893),
                        convert(T, 0.4671734865470852),
                        convert(T, -0.04939742152466368),
                        convert(T, 0.028815162556053812),
                        convert(T, 0.028092718568909044),
                        convert(T, 0.1457093253968254),
                        convert(T, 0.1529431216931217),
                        convert(T, 0.09151785714285714),
                        convert(T, 0.06554705215419501),
                        convert(T, 0.01618992504409171),
                        convert(T, 0.028092718568909044),
                        convert(T, 0.1457093253968254),
                        convert(T, 0.1529431216931217),
                        convert(T, 0.09151785714285714),
                        convert(T, 0.06554705215419501),
                        convert(T, 0.01618992504409171),
                        convert(T, 0.028092718568909044),
                        convert(T, 0.1618992504409171),
                        convert(T, 0.2184901738473167),
                        convert(T, 0.18303571428571427),
                        convert(T, 0.2184901738473167),
                        convert(T, 0.1618992504409171),
                        convert(T, 0.028092718568909044),
                        convert(T, -0.044557986852984274),
                        convert(T, 0.10965442077101599),
                        convert(T, -0.14620589436135464),
                        convert(T, 0.1315853049252192),
                        convert(T, -0.06265966901200913),
                        convert(T, 0.012183824530112887),
                        convert(T, -0.044557986852984274),
                        convert(T, 0.12183824530112887),
                        convert(T, -0.2088655633733638),
                        convert(T, 0.2631706098504384),
                        convert(T, -0.2088655633733638),
                        convert(T, 0.12183824530112887),
                        convert(T, 0.10544201314701572),
                        convert(T, -0.15))
end

struct DPRKN12ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2
    c7::T2
    c8::T2
    c9::T2
    c10::T2
    c11::T2
    c12::T2
    c13::T2
    c14::T2
    c15::T2
    c16::T2
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    a51::T
    # a52::T
    a53::T
    a54::T
    a61::T
    # a62::T
    a63::T
    a64::T
    a65::T
    a71::T
    # a72::T
    a73::T
    a74::T
    a75::T
    a76::T
    a81::T
    # a82::T
    # a83::T
    a84::T
    a85::T
    a86::T
    a87::T
    a91::T
    # a92::T
    a93::T
    a94::T
    a95::T
    a96::T
    a97::T
    a98::T
    a101::T
    # a102::T
    a103::T
    a104::T
    a105::T
    a106::T
    a107::T
    a108::T
    a109::T
    a111::T
    # a112::T
    a113::T
    a114::T
    a115::T
    a116::T
    a117::T
    a118::T
    a119::T
    a1110::T
    a121::T
    # a122::T
    a123::T
    a124::T
    a125::T
    a126::T
    a127::T
    a128::T
    a129::T
    a1210::T
    a1211::T
    a131::T
    # a132::T
    a133::T
    a134::T
    a135::T
    a136::T
    a137::T
    a138::T
    a139::T
    a1310::T
    a1311::T
    a1312::T
    a141::T
    # a142::T
    a143::T
    a144::T
    a145::T
    a146::T
    a147::T
    a148::T
    a149::T
    a1410::T
    a1411::T
    a1412::T
    a1413::T
    a151::T
    # a152::T
    a153::T
    a154::T
    a155::T
    a156::T
    a157::T
    a158::T
    a159::T
    a1510::T
    a1511::T
    a1512::T
    a1513::T
    a1514::T
    a161::T
    # a162::T
    a163::T
    a164::T
    a165::T
    a166::T
    a167::T
    a168::T
    a169::T
    a1610::T
    a1611::T
    a1612::T
    a1613::T
    a1614::T
    a1615::T
    a171::T
    # a172::T
    a173::T
    a174::T
    a175::T
    a176::T
    a177::T
    a178::T
    a179::T
    a1710::T
    a1711::T
    a1712::T
    a1713::T
    a1714::T
    a1715::T
    # a1716::T
    b1::T
    # b2::T
    # b3::T
    # b4::T
    # b5::T
    # b6::T
    b7::T
    b8::T
    b9::T
    b10::T
    b11::T
    b12::T
    b13::T
    b14::T
    b15::T
    # b16::T
    # b17::T
    bp1::T
    # bp2::T
    # bp3::T
    # bp4::T
    # bp5::T
    # bp6::T
    bp7::T
    bp8::T
    bp9::T
    bp10::T
    bp11::T
    bp12::T
    bp13::T
    bp14::T
    bp15::T
    bp16::T
    bp17::T
    btilde1::T
    # btilde2::T
    # btilde3::T
    # btilde4::T
    # btilde5::T
    # btilde6::T
    btilde7::T
    btilde8::T
    btilde9::T
    btilde10::T
    btilde11::T
    btilde12::T
    btilde13::T
    btilde14::T
    btilde15::T
    # btilde16::T
    # btilde17::T
    bptilde1::T
    # bptilde2::T
    # bptilde3::T
    # bptilde4::T
    # bptilde5::T
    # bptilde6::T
    bptilde7::T
    bptilde8::T
    bptilde9::T
    bptilde10::T
    bptilde11::T
    bptilde12::T
    bptilde13::T
    bptilde14::T
    bptilde15::T
    bptilde16::T
    bptilde17::T
end

function DPRKN12ConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 1 // 50)
    c2 = convert(T2, 1 // 25)
    c3 = convert(T2, 1 // 10)
    c4 = convert(T2, 2 // 15)
    c5 = convert(T2, 4 // 25)
    c6 = convert(T2, 1 // 20)
    c7 = convert(T2, 1 // 5)
    c8 = convert(T2, 1 // 4)
    c9 = convert(T2, 1 // 3)
    c10 = convert(T2, 1 // 2)
    c11 = convert(T2, 5 // 9)
    c12 = convert(T2, 3 // 4)
    c13 = convert(T2, 6 // 7)
    c14 = convert(T2, 8437 // 8926)
    c15 = convert(T2, 1)
    c16 = convert(T2, 1)
    a21 = convert(T, 1 // 5000)
    a31 = convert(T, 1 // 3750)
    a32 = convert(T, 1 // 1875)
    a41 = convert(T, 7 // 2400)
    a42 = convert(T, -1 // 240)
    a43 = convert(T, 1 // 160)
    a51 = convert(T, 2 // 1215)
    # a52 = convert(T,0)
    a53 = convert(T, 4 // 729)
    a54 = convert(T, 32 // 18225)
    a61 = convert(T, 152 // 78125)
    # a62 = convert(T,0)
    a63 = convert(T, 1408 // 196875)
    a64 = convert(T, 2048 // 703125)
    a65 = convert(T, 432 // 546875)
    a71 = convert(T, 29 // 51200)
    # a72 = convert(T,0)
    a73 = convert(T, 341 // 387072)
    a74 = convert(T, -151 // 345600)
    a75 = convert(T, 243 // 716800)
    a76 = convert(T, -11 // 110592)
    a81 = convert(T, 37 // 12000)
    # a82 = convert(T,0)
    # a83 = convert(T,0)
    a84 = convert(T, 2 // 1125)
    a85 = convert(T, 27 // 10000)
    a86 = convert(T, 5 // 3168)
    a87 = convert(T, 224 // 20625)
    a91 = convert(T, 100467472123373 // 27511470744477696)
    # a92 = convert(T,0)
    a93 = convert(T, 101066550784375 // 25488568483854336)
    a94 = convert(T, 49478218404275 // 15475202293768704)
    a95 = convert(T, 21990175014231 // 2674726322379776)
    a96 = convert(T, -3576386017671875 // 2723635603703291904)
    a97 = convert(T, 16163228153 // 1654104722787)
    a98 = convert(T, 38747524076705 // 10316801529179136)
    a101 = convert(T, 62178936641284701329 // 16772293867250014666848)
    # a102 = convert(T,0)
    a103 = convert(T, 46108564356250 // 9072835168325103)
    a104 = convert(T, 1522561724950 // 1296119309760729)
    a105 = convert(T, -45978886013453735443 // 2174186242050927827184)
    a106 = convert(T, 299403512366617849203125 // 4981371278573254356053856)
    a107 = convert(T, 15571226634087127616 // 774466927638876610083)
    a108 = convert(T, -133736375367792139885 // 4717207650164066625051)
    a109 = convert(T, 7461389216 // 501451974639)
    a111 = convert(T, 501256914705531962342417557181 // 14270506505142656332600844507392)
    # a112 = convert(T,0)
    a113 = convert(T, -1143766215625 // 132752960853408)
    a114 = convert(T, -6864570325 // 1185294293334)
    a115 = convert(T, 194348369382310456605879163404183 // 99893545535998594328205911551744)
    a116 = convert(T,
                   -94634958447010580589908066176109375 //
                   27549212808177898050085930321520256)
    a117 = convert(T, -17006472665356285286219618514 // 155584463413110817059022733377)
    a118 = convert(T, 33530528814694461893884349656345 // 14270506505142656332600844507392)
    a119 = convert(T, -13439782155791134368 // 17777268379678341919)
    a1110 = convert(T, 1441341768767571 // 13159456712985856)
    a121 = convert(T,
                   parse(BigInt, "105854110734231079069010159870911189747853") //
                   parse(BigInt, "5156624149476760916008179453333467046288864"))
    # a122 = convert(T,0)
    a123 = convert(T, -144579793509250000 // 19842290513127000261)
    a124 = convert(T, -101935644099967250 // 48188419817594143491)
    a125 = convert(T,
                   parse(BigInt, "1585474394319811696785932424388196965") //
                   parse(BigInt, "1709257457318830856936350991091849456"))
    a126 = convert(T,
                   parse(BigInt, "-843499776333774172853009613469456309715703125") //
                   parse(BigInt, "510505790798199330684809765880013237582597536"))
    a127 = convert(T,
                   parse(BigInt, "-15057703799298260121553794369056896088480") //
                   parse(BigInt, "714327132646734138085088291809720015274157"))
    a128 = convert(T,
                   parse(BigInt, "1749840442221344572962864758990584360232600") //
                   parse(BigInt, "1450300542040339007627300471250037606768743"))
    a129 = convert(T, -11255775246405733991656178432768 // 27206626483067760480757659602193)
    a1210 = convert(T, 669010348769579696 // 7368057640845834597)
    a1211 = convert(T, 4598083098752 // 858563707934367)
    a131 = convert(T,
                   parse(BigInt, "-1639758773684715326849438048667467886824967397") //
                   parse(BigInt, "11447568726280607813664651120965112496134881280"))
    # a132 = convert(T,0)
    a133 = convert(T, 3942453384375 // 314673684985856)
    a134 = convert(T, 11737114158175 // 1719466921529856)
    a135 = convert(T,
                   -23710715033675876683332701739887457 //
                   4940189888325748664958546898558976)
    a136 = convert(T,
                   parse(BigInt, "498150575499633273684774666731162498301909124515625") //
                   parse(BigInt, "87415924307623977386706008889913792042985180430336"))
    a137 = convert(T,
                   parse(BigInt, "64881557768202140428371179540010005713998551") //
                   parse(BigInt, "85896810580242200654071863296887242202224768"))
    a138 = convert(T,
                   parse(BigInt, "-2336309182318568698279006266321563486172654055") //
                   parse(BigInt, "18316109962048972501863441793544179993815810048"))
    a139 = convert(T,
                   -493399374030747471036018890494175 // 251658285736841065236836942273664)
    a1310 = convert(T, 418285003077108927126515545155 // 455369916679568501838710898688)
    a1311 = convert(T, -15171723902781457 // 63532954684873728)
    a1312 = convert(T, 1501203688494867 // 9434957026426880)
    a141 = convert(T,
                   parse(BigInt, "34188549803371802849576690267872548602326398788953") //
                   parse(BigInt, "42496542183406636759747616530102745233754251202880"))
    # a142 = convert(T,0)
    a143 = convert(T, -18971246281693750 // 1138830954584356089)
    a144 = convert(T, -59230464334542700 // 2765732318276293359)
    a145 = convert(T,
                   parse(BigInt, "5147939981309774383134903239728881770043") //
                   parse(BigInt, "305929030949718561059100251282184099064"))
    a146 = convert(T,
                   parse(BigInt,
                         "-3625720213550267723370658302114678215563058405229078120") //
                   parse(BigInt, "324512095420929759624784749347170583153994213035432256"))
    a147 = convert(T,
                   parse(BigInt, "-60305503318319653518547439098565661266182518307816") //
                   parse(BigInt, "17856872599361492097414471889911176856851308259643"))
    a148 = convert(T,
                   parse(BigInt, "-1036461878759982363277481306266144563833492657780645") //
                   parse(BigInt, "67994467493450618815596186448164392374006801924608"))
    a149 = convert(T,
                   parse(BigInt, "128398681100219349205889126776607047000") //
                   parse(BigInt, "7473801441221286756994805323613917077"))
    a1410 = convert(T, -49156374556350058671822606102117 // 9039888303968618912866414995904)
    a1411 = convert(T, 12253036339964386945 // 8828680926314891943)
    a1412 = convert(T, -647188390508758231059 // 1092148506009694282240)
    a1413 = convert(T, 10915833599872 // 368729913707897)
    a151 = convert(T,
                   parse(BigInt,
                         "-4939337286263213195547765488387521892799075623007291241961609516532") //
                   parse(BigInt,
                         "5408250052307451520718178852915698257207815452080611897685945761264"))
    # a152 = convert(T,0)
    a153 = convert(T,
                   7588799849596321243074032368290625 //
                   parse(BigInt, "3147217749590114939838670370597819616"))
    a154 = convert(T,
                   16870665568420512953501332587233725 //
                   955405388268427749593882076788623812)
    a155 = convert(T,
                   parse(BigInt,
                         "-808642515918378014850308582271476014669568437579087796060") //
                   parse(BigInt,
                         "54447992506702009927986632715967769032585338753056786562"))
    a156 = convert(T,
                   parse(BigInt,
                         "4610328329649866588704236006423149172472141907645890762410296050212") //
                   parse(BigInt,
                         "2135428689710103309390449198881479603148467934048051598947383737508"))
    a157 = convert(T,
                   parse(BigInt,
                         "4159963831215576225909381034291748993887819834160487158570788681") //
                   parse(BigInt,
                         "1040533184037697645660563795162185415624171583014576682740416336"))
    a158 = convert(T,
                   parse(BigInt,
                         "7381392142124351279433801934148706553542137071890521365664606664449580") //
                   parse(BigInt,
                         "259596002510757672994472584939953516345975141699869371088925396540699"))
    a159 = convert(T,
                   parse(BigInt,
                         "-3336834334584052813468828675971359774694437229547862706920") //
                   parse(BigInt,
                         "132102862435303266640535426836147775872819092781208127980"))
    a1510 = convert(T,
                    parse(BigInt,
                          "426619379967412086875039012957475466130081426048213491790") //
                    parse(BigInt,
                          "55162410119399855550108207148248549410926885937244965785"))
    a1511 = convert(T,
                    parse(BigInt, "-630755628691078947314733435975762542732598947") //
                    parse(BigInt, "333503232300511886435069380727586592765317456"))
    a1512 = convert(T,
                    parse(BigInt, "1522350657470125698997653827133798314909646891") //
                    parse(BigInt, "1520094067152619944607524353149267399623188480"))
    a1513 = convert(T,
                    305575414262755427083262606101825880 //
                    parse(BigInt, "65839748482572312891297405431209259829"))
    a1514 = convert(T,
                    parse(BigInt, "256624643108055110568255672032710477795") //
                    parse(BigInt, "22874609758516552135947898572671559986304"))
    a161 = convert(T,
                   parse(BigInt,
                         "-571597862947184314270186718640978947715678864684269066846") //
                   parse(BigInt,
                         "2077055064880303907616135969012720011907767004397744786340"))
    # a162 = convert(T,0)
    a163 = convert(T, 66981514290625 // 1829501741761029)
    a164 = convert(T, 43495576635800 // 4443075658562499)
    a165 = convert(T,
                   -127865248353371207265315478623656127 //
                   10401415428935853634424440540325344)
    a166 = convert(T,
                   parse(BigInt,
                         "1316565142658075739557231574080234814338066993483960326560") //
                   parse(BigInt,
                         "92668695535091962564795912774190176478892159517481612467"))
    a167 = convert(T,
                   parse(BigInt,
                         "3881494143728609118531066904799685950051960514138645179820") //
                   parse(BigInt,
                         "2446349095978358868919950548516272963929118212742344026549"))
    a168 = convert(T,
                   parse(BigInt,
                         "162922667049680755852592453758428194006198229544701786842910") //
                   parse(BigInt,
                         "66288722243155885736983218667976563740242178853010092663614"))
    a169 = convert(T,
                   parse(BigInt, "-43986024977384568043684084266385512680544563954") //
                   parse(BigInt, "4922783599524658241955780540171948284522386185"))
    a1610 = convert(T,
                    parse(BigInt, "285912200202585226675651763671663063668290787") //
                    parse(BigInt, "65371192072964016939690070594254881767827200"))
    a1611 = convert(T, -6776815256667778089672518929 // 3693654613173093729492918708)
    a1612 = convert(T,
                    398946554885847045598775476868169 // 344154261237450078839899047372800)
    a1613 = convert(T, -76630698033396272 // 4432017119727044925)
    a1614 = convert(T, 28401702316003037 // 1469612686944417840)
    a1615 = convert(T,
                    66049942462586341419969330578128801 //
                    parse(BigInt, "12691068622536592094919763114637498325"))
    a171 = convert(T,
                   parse(BigInt,
                         "83940754497395557520874219603241359529066454343054832302344735") //
                   parse(BigInt,
                         "64192596456995578553872477759926464976144474354415663868673233"))
    # a172 = convert(T,0)
    a173 = convert(T, 892543892035485503125 // 51401651664490002607536)
    a174 = convert(T, -12732238157949399705325 // 686579204375687891972088)
    a175 = convert(T,
                   parse(BigInt, "5290376174838819557032232941734928484252549") //
                   parse(BigInt, "357179779572898187570048915214361602000384"))
    a176 = convert(T,
                   parse(BigInt,
                         "26873229338017506937199991804717456666650215387938173031932210") //
                   parse(BigInt,
                         "2863980005760296740624015421425947092438943496681472214589916"))
    a177 = convert(T,
                   parse(BigInt,
                         "-1976497866818803305857417297961598735637414137241493515492778650") //
                   parse(BigInt,
                         "378029217824623393200881653405474359138017953416246216408422692"))
    a178 = convert(T,
                   parse(BigInt,
                         "-1002860756304839757040188283199900676042073362417943601440986856950") //
                   parse(BigInt,
                         "20486915674765670626893195919603679319429068544972409068469849579"))
    a179 = convert(T,
                   parse(BigInt,
                         "87398661196965758104117684348440686081062878816711392590") //
                   parse(BigInt, "2282122412587168891929052689609009868137678763277087160"))
    a1710 = convert(T,
                    parse(BigInt,
                          "-7922242431969626895355493632206885458496418610471389") //
                    parse(BigInt, "748272134517487495468365669337985635214015258726400"))
    a1711 = convert(T,
                    parse(BigInt, "2777643183645212014464950387658055285") //
                    parse(BigInt, "1141545470045611737197667093465955392"))
    a1712 = convert(T,
                    parse(BigInt, "-1372659703515496442825084239977218110461") //
                    parse(BigInt, "1313121960368535725613950174847107891200"))
    a1713 = convert(T, 6144417902699179309851023 // 85608793932459282773805825)
    a1714 = convert(T, 140294243355138853053241 // 64884622846351585391642880)
    a1715 = convert(T,
                    parse(BigInt, "168671028523891369934964082754523881107337") //
                    parse(BigInt, "24062875279623260368388427013982199424119600"))
    # a1716 = convert(T,0)
    b1 = convert(T, 63818747 // 5262156900)
    # b2 = convert(T,0)
    # b3 = convert(T,0)
    # b4 = convert(T,0)
    # b5 = convert(T,0)
    # b6 = convert(T,0)
    b7 = convert(T, 22555300000000 // 261366897038247)
    b8 = convert(T, 1696514453125 // 6717619827072)
    b9 = convert(T, -45359872 // 229764843)
    b10 = convert(T, 19174962087 // 94371046000)
    b11 = convert(T, -19310468 // 929468925)
    b12 = convert(T, 16089185487681 // 146694672924800)
    b13 = convert(T, 1592709632 // 41841694125)
    b14 = convert(T, 52675701958271 // 4527711056573100)
    b15 = convert(T,
                  parse(BigInt, "12540904472870916741199505796420811396") //
                  parse(BigInt, "2692319557780977037279406889319526430375"))
    # b16 = convert(T,0)
    # b17 = convert(T,0)
    bp1 = convert(T, 63818747 // 5262156900)
    # bp2 = convert(T,0)
    # bp3 = convert(T,0)
    # bp4 = convert(T,0)
    # bp5 = convert(T,0)
    # bp6 = convert(T,0)
    bp7 = convert(T, 451106000000000 // 4965971043726693)
    bp8 = convert(T, 8482572265625 // 26870479308288)
    bp9 = convert(T, -181439488 // 689294529)
    bp10 = convert(T, 57524886261 // 188742092000)
    bp11 = convert(T, -38620936 // 929468925)
    bp12 = convert(T, 144802669389129 // 586778691699200)
    bp13 = convert(T, 6370838528 // 41841694125)
    bp14 = convert(T, 368729913707897 // 4527711056573100)
    bp15 = convert(T,
                   parse(BigInt, "111940113324845802831946788738852162520696") //
                   parse(BigInt, "1316544263754897771229629968877248424453375"))
    bp16 = convert(T, -113178587 // 12362232960)
    bp17 = convert(T, 1 // 40)

    btilde1 = convert(T,
                      Int64(63818747) // Int64(5262156900) -
                      Int64(27121957) // Int64(1594593000))
    # btilde2 = convert(T,0)
    # btilde3 = convert(T,0)
    # btilde4 = convert(T,0)
    # btilde5 = convert(T,0)
    # btilde6 = convert(T,0)
    btilde7 = convert(T,
                      Int64(22555300000000) // Int64(261366897038247) -
                      Int64(4006163300000) // Int64(55441463008113))
    btilde8 = convert(T,
                      Int64(1696514453125) // Int64(6717619827072) -
                      Int64(9466403125) // Int64(25445529648))
    btilde9 = convert(T,
                      Int64(-45359872) // Int64(229764843) +
                      Int64(163199648) // Int64(406149975))
    btilde10 = convert(T,
                       Int64(19174962087) // Int64(94371046000) -
                       Int64(23359833) // Int64(69636250))
    btilde11 = convert(T,
                       Int64(-19310468) // Int64(929468925) +
                       Int64(18491714) // Int64(140828625))
    btilde12 = convert(T,
                       Int64(16089185487681) // Int64(146694672924800) -
                       Int64(11052304606701) // Int64(58344472186000))
    btilde13 = convert(T,
                       Int64(1592709632) // Int64(41841694125) -
                       Int64(1191129152) // Int64(44377554375))
    btilde14 = convert(T,
                       Int64(52675701958271) // Int64(4527711056573100) -
                       Int64(2033811086741) // Int64(124730332137000))
    btilde15 = convert(T,
                       parse(BigInt, "12540904472870916741199505796420811396") //
                       parse(BigInt, "2692319557780977037279406889319526430375") -
                       parse(BigInt, "3616943474975740389660406409450169802") //
                       parse(BigInt, "951830146690244407118982233597812374375"))
    # btilde16 = convert(T,0)
    # btilde17 = convert(T,0)
    bptilde1 = convert(T,
                       Int64(63818747) // Int64(5262156900) -
                       Int64(27121957) // Int64(1594593000))
    # bptilde2 = convert(T,0)
    # bptilde3 = convert(T,0)
    # bptilde4 = convert(T,0)
    # bptilde5 = convert(T,0)
    # bptilde6 = convert(T,0)
    bptilde7 = convert(T,
                       Int64(451106000000000) // Int64(4965971043726693) -
                       Int64(4217014000000) // Int64(55441463008113))
    bptilde8 = convert(T,
                       Int64(8482572265625) // Int64(26870479308288) -
                       Int64(47332015625) // Int64(101782118592))
    bptilde9 = convert(T,
                       Int64(-181439488) // Int64(689294529) +
                       Int64(652798592) // Int64(1218449925))
    bptilde10 = convert(T,
                        Int64(57524886261) // Int64(188742092000) -
                        Int64(70079499) // Int64(139272500))
    bptilde11 = convert(T,
                        Int64(-38620936) // Int64(929468925) +
                        Int64(36983428) // Int64(140828625))
    bptilde12 = convert(T,
                        Int64(144802669389129) // Int64(586778691699200) -
                        Int64(99470741460309) // Int64(233377888744000))
    bptilde13 = convert(T,
                        Int64(6370838528) // Int64(41841694125) -
                        Int64(4764516608) // Int64(44377554375))
    bptilde14 = convert(T,
                        Int64(368729913707897) // Int64(4527711056573100) -
                        Int64(14236677607187) // Int64(124730332137000))
    bptilde15 = convert(T,
                        parse(BigInt, "111940113324845802831946788738852162520696") //
                        parse(BigInt, "1316544263754897771229629968877248424453375") -
                        parse(BigInt, "198066487470143918516004831967805004004") //
                        parse(BigInt, "2855490440070733221356946700793437123125"))
    bptilde16 = convert(T, Int64(-113178587) // Int64(12362232960) - Int64(1) // Int64(50))
    bptilde17 = convert(T, 1 // 40)
    DPRKN12ConstantCache(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15,
                         c16, a21, a31, a32, a41, a42, a43, a51, a53, a54, a61, a63, a64,
                         a65, a71, a73, a74, a75, a76, a81, a84, a85, a86, a87, a91, a93,
                         a94, a95, a96, a97, a98, a101, a103, a104, a105, a106, a107, a108,
                         a109, a111, a113, a114, a115, a116, a117, a118, a119, a1110, a121,
                         a123, a124, a125, a126, a127, a128, a129, a1210, a1211, a131, a133,
                         a134, a135, a136, a137, a138, a139, a1310, a1311, a1312, a141,
                         a143, a144, a145, a146, a147, a148, a149, a1410, a1411, a1412,
                         a1413, a151, a153, a154, a155, a156, a157, a158, a159, a1510,
                         a1511, a1512, a1513, a1514, a161, a163, a164, a165, a166, a167,
                         a168, a169, a1610, a1611, a1612, a1613, a1614, a1615, a171, a173,
                         a174, a175, a176, a177, a178, a179, a1710, a1711, a1712, a1713,
                         a1714, a1715, b1, b7, b8, b9, b10, b11, b12, b13, b14, b15, bp1,
                         bp7, bp8, bp9, bp10, bp11, bp12, bp13, bp14, bp15, bp16, bp17,
                         btilde1, btilde7, btilde8, btilde9, btilde10, btilde11, btilde12,
                         btilde13, btilde14, btilde15, bptilde1, bptilde7, bptilde8,
                         bptilde9, bptilde10, bptilde11, bptilde12, bptilde13, bptilde14,
                         bptilde15, bptilde16, bptilde17)
end

function DPRKN12ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    DPRKN12ConstantCache(convert(T2, 2.0e-2),
                         convert(T2, 4.0e-2),
                         convert(T2, 1.0e-1),
                         convert(T2, 1.33333333333333333333333333333e-1),
                         convert(T2, 1.6e-1),
                         convert(T2, 5.0e-2),
                         convert(T2, 2.0e-1),
                         convert(T2, 2.5e-1),
                         convert(T2, 3.33333333333333333333333333333e-1),
                         convert(T2, 5.0e-1),
                         convert(T2, 5.55555555555555555555555555556e-1),
                         convert(T2, 7.5e-1),
                         convert(T2, 8.57142857142857142857142857143e-1),
                         convert(T2, 9.45216222272014340129957427739e-1),
                         convert(T2, 1.0e0),
                         convert(T2, 1.0e0),
                         convert(T, 2.0e-4),
                         convert(T, 2.66666666666666666666666666667e-4),
                         convert(T, 5.33333333333333333333333333333e-4),
                         convert(T, 2.91666666666666666666666666667e-3),
                         convert(T, -4.16666666666666666666666666667e-3),
                         convert(T, 6.25e-3),
                         convert(T, 1.64609053497942386831275720165e-3),
                         convert(T, 5.48696844993141289437585733882e-3),
                         convert(T, 1.75582990397805212620027434842e-3),
                         convert(T, 1.9456e-3),
                         convert(T, 7.15174603174603174603174603175e-3),
                         convert(T, 2.91271111111111111111111111111e-3),
                         convert(T, 7.89942857142857142857142857143e-4),
                         convert(T, 5.6640625e-4),
                         convert(T, 8.80973048941798941798941798942e-4),
                         convert(T, -4.36921296296296296296296296296e-4),
                         convert(T, 3.39006696428571428571428571429e-4),
                         convert(T, -9.94646990740740740740740740741e-5),
                         convert(T, 3.08333333333333333333333333333e-3),
                         convert(T, 1.77777777777777777777777777778e-3),
                         convert(T, 2.7e-3),
                         convert(T, 1.57828282828282828282828282828e-3),
                         convert(T, 1.08606060606060606060606060606e-2),
                         convert(T, 3.65183937480112971375119150338e-3),
                         convert(T, 3.96517171407234306617557289807e-3),
                         convert(T, 3.19725826293062822350093426091e-3),
                         convert(T, 8.22146730685543536968701883401e-3),
                         convert(T, -1.31309269595723798362013884863e-3),
                         convert(T, 9.77158696806486781562609494147e-3),
                         convert(T, 3.75576906923283379487932641079e-3),
                         convert(T, 3.70724106871850081019565530521e-3),
                         convert(T, 5.08204585455528598076108163479e-3),
                         convert(T, 1.17470800217541204473569104943e-3),
                         convert(T, -2.11476299151269914996229766362e-2),
                         convert(T, 6.01046369810788081222573525136e-2),
                         convert(T, 2.01057347685061881846748708777e-2),
                         convert(T, -2.83507501229335808430366774368e-2),
                         convert(T, 1.48795689185819327555905582479e-2),
                         convert(T, 3.51253765607334415311308293052e-2),
                         convert(T, -8.61574919513847910340576078545e-3),
                         convert(T, -5.79144805100791652167632252471e-3),
                         convert(T, 1.94555482378261584239438810411e0),
                         convert(T, -3.43512386745651359636787167574e0),
                         convert(T, -1.09307011074752217583892572001e-1),
                         convert(T, 2.3496383118995166394320161088e0),
                         convert(T, -7.56009408687022978027190729778e-1),
                         convert(T, 1.09528972221569264246502018618e-1),
                         convert(T, 2.05277925374824966509720571672e-2),
                         convert(T, -7.28644676448017991778247943149e-3),
                         convert(T, -2.11535560796184024069259562549e-3),
                         convert(T, 9.27580796872352224256768033235e-1),
                         convert(T, -1.65228248442573667907302673325e0),
                         convert(T, -2.10795630056865698191914366913e-2),
                         convert(T, 1.20653643262078715447708832536e0),
                         convert(T, -4.13714477001066141324662463645e-1),
                         convert(T, 9.07987398280965375956795739516e-2),
                         convert(T, 5.35555260053398504916870658215e-3),
                         convert(T, -1.43240788755455150458921091632e-1),
                         convert(T, 1.25287037730918172778464480231e-2),
                         convert(T, 6.82601916396982712868112411737e-3),
                         convert(T, -4.79955539557438726550216254291e0),
                         convert(T, 5.69862504395194143379169794156e0),
                         convert(T, 7.55343036952364522249444028716e-1),
                         convert(T, -1.27554878582810837175400796542e-1),
                         convert(T, -1.96059260511173843289133255423e0),
                         convert(T, 9.18560905663526240976234285341e-1),
                         convert(T, -2.38800855052844310534827013402e-1),
                         convert(T, 1.59110813572342155138740170963e-1),
                         convert(T, 8.04501920552048948697230778134e-1),
                         convert(T, -1.66585270670112451778516268261e-2),
                         convert(T, -2.1415834042629734811731437191e-2),
                         convert(T, 1.68272359289624658702009353564e1),
                         convert(T, -1.11728353571760979267882984241e1),
                         convert(T, -3.37715929722632374148856475521e0),
                         convert(T, -1.52433266553608456461817682939e1),
                         convert(T, 1.71798357382154165620247684026e1),
                         convert(T, -5.43771923982399464535413738556e0),
                         convert(T, 1.38786716183646557551256778839e0),
                         convert(T, -5.92582773265281165347677029181e-1),
                         convert(T, 2.96038731712973527961592794552e-2),
                         convert(T, -9.13296766697358082096250482648e-1),
                         convert(T, 2.41127257578051783924489946102e-3),
                         convert(T, 1.76581226938617419820698839226e-2),
                         convert(T, -1.48516497797203838246128557088e1),
                         convert(T, 2.15897086700457560030782161561e0),
                         convert(T, 3.99791558311787990115282754337e0),
                         convert(T, 2.84341518002322318984542514988e1),
                         convert(T, -2.52593643549415984378843352235e1),
                         convert(T, 7.7338785423622373655340014114e0),
                         convert(T, -1.8913028948478674610382580129e0),
                         convert(T, 1.00148450702247178036685959248e0),
                         convert(T, 4.64119959910905190510518247052e-3),
                         convert(T, 1.12187550221489570339750499063e-2),
                         convert(T, -2.75196297205593938206065227039e-1),
                         convert(T, 3.66118887791549201342293285553e-2),
                         convert(T, 9.7895196882315626246509967162e-3),
                         convert(T, -1.2293062345886210304214726509e1),
                         convert(T, 1.42072264539379026942929665966e1),
                         convert(T, 1.58664769067895368322481964272e0),
                         convert(T, 2.45777353275959454390324346975e0),
                         convert(T, -8.93519369440327190552259086374e0),
                         convert(T, 4.37367273161340694839327077512e0),
                         convert(T, -1.83471817654494916304344410264e0),
                         convert(T, 1.15920852890614912078083198373e0),
                         convert(T, -1.72902531653839221518003422953e-2),
                         convert(T, 1.93259779044607666727649875324e-2),
                         convert(T, 5.20444293755499311184926401526e-3),
                         convert(T, 1.30763918474040575879994562983e0),
                         convert(T, 1.73641091897458418670879991296e-2),
                         convert(T, -1.8544456454265795024362115588e-2),
                         convert(T, 1.48115220328677268968478356223e1),
                         convert(T, 9.38317630848247090787922177126e0),
                         convert(T, -5.2284261999445422541474024553e0),
                         convert(T, -4.89512805258476508040093482743e1),
                         convert(T, 3.82970960343379225625583875836e1),
                         convert(T, -1.05873813369759797091619037505e1),
                         convert(T, 2.43323043762262763585119618787e0),
                         convert(T, -1.04534060425754442848652456513e0),
                         convert(T, 7.17732095086725945198184857508e-2),
                         convert(T, 2.16221097080827826905505320027e-3),
                         convert(T, 7.00959575960251423699282781988e-3),
                         convert(T, 0.012127868517185414),
                         convert(T, 0.08629746251568875),
                         convert(T, 0.2525469581187147),
                         convert(T, -0.1974186799326823),
                         convert(T, 0.2031869190789726),
                         convert(T, -0.020775808077714918),
                         convert(T, 0.10967804874502014),
                         convert(T, 0.038065132526466504),
                         convert(T, 0.01163406880432423),
                         convert(T, 0.0046580297040248785),
                         convert(T, 0.012127868517185414),
                         convert(T, 0.09083943422704079),
                         convert(T, 0.3156836976483934),
                         convert(T, -0.2632249065769097),
                         convert(T, 0.3047803786184589),
                         convert(T, -0.041551616155429835),
                         convert(T, 0.2467756096762953),
                         convert(T, 0.15226053010586602),
                         convert(T, 0.08143848163026961),
                         convert(T, 0.08502571193890811),
                         convert(T, -0.009155189630077963),
                         convert(T, 0.025),
                         convert(T, -0.004880833389821577),
                         convert(T, 0.014038126584857338),
                         convert(T, -0.11947921920803832),
                         convert(T, 0.20440246507662121),
                         convert(T, -0.1322681492223791),
                         convert(T, 0.11053069299761689),
                         convert(T, -0.07975385787102851),
                         convert(T, 0.011224330486437457),
                         convert(T, -0.004671596801593694),
                         convert(T, 0.0008580413473282841),
                         convert(T, -0.004880833389821577),
                         convert(T, 0.014776975352481408),
                         convert(T, -0.1493490240100479),
                         convert(T, 0.2725366201021616),
                         convert(T, -0.19840222383356862),
                         convert(T, 0.22106138599523378),
                         convert(T, -0.17944618020981415),
                         convert(T, 0.04489732194574983),
                         convert(T, -0.03270117761115586),
                         convert(T, 0.015662325288859434),
                         convert(T, -0.029155189630077964),
                         convert(T, 0.025))
end
