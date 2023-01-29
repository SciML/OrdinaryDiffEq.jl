struct BS3ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    a21::T
    a32::T
    a41::T
    a42::T
    a43::T
    c1::T2
    c2::T2
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
end

"""
constructBogakiShampine3()

Constructs the tableau object for the Bogakai-Shampine Order 2/3 method.
"""
function BS3ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    a21 = convert(T, 0.5)
    a32 = convert(T, 0.75)
    a41 = convert(T, 0.2222222222222222)
    a42 = convert(T, 0.3333333333333333)
    a43 = convert(T, 0.4444444444444444)
    c1 = convert(T2, 0.5)
    c2 = convert(T2, 0.75)
    # b1 = convert(T,0.2916666666666667)
    # b2 = convert(T,0.25)
    # b3 = convert(T,0.3333333333333333)
    # b4 = convert(T,0.125)
    btilde1 = convert(T, 0.06944444444444445)
    btilde2 = convert(T, -0.08333333333333333)
    btilde3 = convert(T, -0.1111111111111111)
    btilde4 = convert(T, 0.125)
    BS3ConstantCache(a21, a32, a41, a42, a43, c1, c2, btilde1, btilde2, btilde3, btilde4)
end

"""
constructBogakiShampine3()

Constructs the tableau object for the Bogakai-Shampine Order 2/3 method.
"""
function BS3ConstantCache(T::Type, T2::Type)
    a21 = convert(T, 1 // 2)
    a32 = convert(T, 3 // 4)
    a41 = convert(T, 2 // 9)
    a42 = convert(T, 1 // 3)
    a43 = convert(T, 4 // 9)
    c1 = convert(T2, 1 // 2)
    c2 = convert(T2, 3 // 4)
    # b1 = convert(T,7//24)
    # b2 = convert(T,1//4)
    # b3 = convert(T,1//3)
    # b4 = convert(T,1//8)
    btilde1 = convert(T, 5 // 72)
    btilde2 = convert(T, -1 // 12)
    btilde3 = convert(T, -1 // 9)
    btilde4 = convert(T, 1 // 8)
    BS3ConstantCache(a21, a32, a41, a42, a43, c1, c2, btilde1, btilde2, btilde3, btilde4)
end

struct OwrenZen3ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    c1::T2
    c2::T2
    btilde1::T
    btilde2::T
    btilde3::T
    r13::T
    r12::T
    r23::T
    r22::T
    r33::T
    r32::T
end

function OwrenZen3ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    a21 = convert(T, 0.5217391304347826)
    a31 = convert(T, -0.18133333333333335)
    a32 = convert(T, 0.9813333333333333)
    a41 = convert(T, 0.2152777777777778)
    a42 = convert(T, 0.4592013888888889)
    a43 = convert(T, 0.3255208333333333)
    c1 = convert(T2, 0.5217391304347826)
    c2 = convert(T2, 0.8)
    # b1 = convert(T,0.041666666666666664)
    # b2 = convert(T,0.9583333333333334)
    btilde1 = convert(T, -0.1736111111111111)
    btilde2 = convert(T, 0.4991319444444444)
    btilde3 = convert(T, -0.3255208333333333)
    r13 = convert(T, 0.5694444444444444)
    r12 = convert(T, -1.3541666666666667)
    r23 = convert(T, -0.9184027777777778)
    r22 = convert(T, 1.3776041666666667)
    r33 = convert(T, -0.6510416666666666)
    r32 = convert(T, 0.9765625)
    OwrenZen3ConstantCache(a21, a31, a32, a41, a42, a43, c1, c2, btilde1, btilde2, btilde3,
                           r13, r12, r23, r22, r33, r32)
end

function OwrenZen3ConstantCache(T, T2)
    a21 = convert(T, 12 // 23)
    a31 = convert(T, -68 // 375)
    a32 = convert(T, 368 // 375)
    a41 = convert(T, 31 // 144)
    a42 = convert(T, 529 // 1152)
    a43 = convert(T, 125 // 384)
    c1 = convert(T2, 12 // 23)
    c2 = convert(T2, 4 // 5)
    # b1 = convert(T,1//24)
    # b2 = convert(T,23//24)
    btilde1 = convert(T, -25 // 144)
    btilde2 = convert(T, 575 // 1152)
    btilde3 = convert(T, -125 // 384)
    r13 = convert(T, 41 // 72)
    r12 = convert(T, -65 // 48)
    r23 = convert(T, -529 // 576)
    r22 = convert(T, 529 // 384)
    r33 = convert(T, -125 // 192)
    r32 = convert(T, 125 // 128)
    OwrenZen3ConstantCache(a21, a31, a32, a41, a42, a43, c1, c2, btilde1, btilde2, btilde3,
                           r13, r12, r23, r22, r33, r32)
end

struct OwrenZen4ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
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
    a63::T
    a64::T
    a65::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    btilde1::T
    btilde3::T
    btilde4::T
    btilde5::T
    r14::T
    r13::T
    r12::T
    r34::T
    r33::T
    r32::T
    r44::T
    r43::T
    r42::T
    r54::T
    r53::T
    r52::T
    r64::T
    r63::T
    r62::T
end

function OwrenZen4ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    a21 = convert(T, 0.16666666666666666)
    a31 = convert(T, 0.03214024835646457)
    a32 = convert(T, 0.26515704894083275)
    a41 = convert(T, 0.6895990230002036)
    a42 = convert(T, -1.6993690209647874)
    a43 = convert(T, 1.6568288214939955)
    a51 = convert(T, -0.09002509947964493)
    a52 = convert(T, 0.6817777777777778)
    a53 = convert(T, -0.2402791551882461)
    a54 = convert(T, 0.5151931435567799)
    a61 = convert(T, 0.08990252172070354)
    a63 = convert(T, 0.4360623278236915)
    a64 = convert(T, 0.1842858372687918)
    a65 = convert(T, 0.2897493131868132)
    c1 = convert(T2, 0.16666666666666666)
    c2 = convert(T2, 0.2972972972972973)
    c3 = convert(T2, 0.6470588235294118)
    c4 = convert(T2, 0.8666666666666667)
    # b1 = convert(T,0.27823691460055094)
    # b3 = convert(T,-0.09428374655647383)
    # b4 = convert(T,0.8160468319559229)
    btilde1 = convert(T, 0.18833439287984743)
    btilde3 = convert(T, -0.5303460743801653)
    btilde4 = convert(T, 0.6317609946871311)
    btilde5 = convert(T, -0.2897493131868132)
    r14 = convert(T, -1.0513495872621479)
    r13 = convert(T, 2.922894131082889)
    r12 = convert(T, -2.7816420221000375)
    r34 = convert(T, 2.4266369235379472)
    r33 = convert(T, -5.725398502723277)
    r32 = convert(T, 3.7348239070090217)
    r44 = convert(T, -0.7705135914137909)
    r43 = convert(T, 1.1724555082899983)
    r42 = convert(T, -0.21765607960741548)
    r54 = convert(T, -2.8643157295948325)
    r53 = convert(T, 5.149132832816039)
    r52 = convert(T, -1.9950677900343932)
    r64 = convert(T, 2.2595419847328246)
    r63 = convert(T, -3.519083969465649)
    r62 = convert(T, 1.2595419847328244)
    OwrenZen4ConstantCache(a21, a31, a32, a41, a42, a43, a51, a52, a53, a54,
                           a61, a63, a64, a65, c1, c2, c3, c4, btilde1, btilde3, btilde4,
                           btilde5,
                           r14, r13, r12, r34, r33, r32, r44, r43, r42,
                           r54, r53, r52, r64, r63, r62)
end

function OwrenZen4ConstantCache(T, T2)
    a21 = convert(T, 1 // 6)
    a31 = convert(T, 44 // 1369)
    a32 = convert(T, 363 // 1369)
    a41 = convert(T, 3388 // 4913)
    a42 = convert(T, -8349 // 4913)
    a43 = convert(T, 8140 // 4913)
    a51 = convert(T, -36764 // 408375)
    a52 = convert(T, 767 // 1125)
    a53 = convert(T, -32708 // 136125)
    a54 = convert(T, 210392 // 408375)
    a61 = convert(T, 1697 // 18876)
    a63 = convert(T, 50653 // 116160)
    a64 = convert(T, 299693 // 1626240)
    a65 = convert(T, 3375 // 11648)
    c1 = convert(T2, 1 // 6)
    c2 = convert(T2, 11 // 37)
    c3 = convert(T2, 11 // 17)
    c4 = convert(T2, 13 // 15)
    # b1 = convert(T,101//363)
    # b3 = convert(T,-1369//14520)
    # b4 = convert(T,11849//14520)
    btilde1 = convert(T, 1185 // 6292)
    btilde3 = convert(T, -4107 // 7744)
    btilde4 = convert(T, 68493 // 108416)
    btilde5 = convert(T, -3375 // 11648)
    r14 = convert(T, -866577 // 824252)
    r13 = convert(T, 1806901 // 618189)
    r12 = convert(T, -104217 // 37466)
    r34 = convert(T, 12308679 // 5072320)
    r33 = convert(T, -2178079 // 380424)
    r32 = convert(T, 861101 // 230560)
    r44 = convert(T, -7816583 // 10144640)
    r43 = convert(T, 6244423 // 5325936)
    r42 = convert(T, -63869 // 293440)
    r54 = convert(T, -624375 // 217984)
    r53 = convert(T, 982125 // 190736)
    r52 = convert(T, -1522125 // 762944)
    r64 = convert(T, 296 // 131)
    r63 = convert(T, -461 // 131)
    r62 = convert(T, 165 // 131)
    OwrenZen4ConstantCache(a21, a31, a32, a41, a42, a43, a51, a52, a53, a54,
                           a61, a63, a64, a65, c1, c2, c3, c4, btilde1, btilde3, btilde4,
                           btilde5,
                           r14, r13, r12, r34, r33, r32, r44, r43, r42,
                           r54, r53, r52, r64, r63, r62)
end

struct OwrenZen5ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    a21::T
    a31::T
    a32::T
    a41::T
    a42::T
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
    a83::T
    a84::T
    a85::T
    a86::T
    a87::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2
    btilde1::T
    btilde3::T
    btilde4::T
    btilde5::T
    btilde6::T
    btilde7::T
    r15::T
    r14::T
    r13::T
    r12::T
    r35::T
    r34::T
    r33::T
    r32::T
    r45::T
    r44::T
    r43::T
    r42::T
    r55::T
    r54::T
    r53::T
    r52::T
    r65::T
    r64::T
    r63::T
    r62::T
    r75::T
    r74::T
    r73::T
    r72::T
    r85::T
    r84::T
    r83::T
    r82::T
end

function OwrenZen5ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    a21 = convert(T, 0.16666666666666666)
    a31 = convert(T, 0.0625)
    a32 = convert(T, 0.1875)
    a41 = convert(T, 0.25)
    a42 = convert(T, -0.75)
    a51 = convert(T, -0.75)
    a52 = convert(T, 3.75)
    a53 = convert(T, -3.0)
    a54 = convert(T, 0.5)
    a61 = convert(T, 0.26895043731778423)
    a62 = convert(T, -0.7084548104956269)
    a63 = convert(T, 0.8658892128279884)
    a64 = convert(T, 0.15462307371928363)
    a65 = convert(T, 0.06184922948771345)
    a71 = convert(T, -0.02947695035460993)
    a72 = convert(T, 0.18500664893617022)
    a73 = convert(T, 0.4802345261121857)
    a74 = convert(T, -0.5337849069148937)
    a75 = convert(T, -0.013090093085106383)
    a76 = convert(T, 0.7861107753062541)
    a81 = convert(T, 0.08783068783068783)
    a83 = convert(T, 0.3006060606060606)
    a84 = convert(T, 0.22777777777777777)
    a85 = convert(T, 0.027777777777777776)
    a86 = convert(T, 0.06218596218596219)
    a87 = convert(T, 0.2938217338217338)
    c1 = convert(T2, 0.16666666666666666)
    c2 = convert(T2, 0.25)
    c3 = convert(T2, 0.5)
    c4 = convert(T2, 0.5)
    c5 = convert(T2, 0.6428571428571429)
    c6 = convert(T2, 0.875)
    # b1 = convert(T,-0.1111111111111111)
    # b3 = convert(T,1.2121212121212122)
    # b4 = convert(T,-1.75)
    # b5 = convert(T,-0.08333333333333333)
    # b6 = convert(T,1.7323232323232323)
    btilde1 = convert(T, -0.19894179894179895)
    btilde3 = convert(T, 0.9115151515151515)
    btilde4 = convert(T, -1.9777777777777779)
    btilde5 = convert(T, -0.1111111111111111)
    btilde6 = convert(T, 1.67013727013727)
    btilde7 = convert(T, -0.2938217338217338)
    r15 = convert(T, 1.892063492063492)
    r14 = convert(T, -6.067155067155067)
    r13 = convert(T, 7.282458282458283)
    r12 = convert(T, -4.0195360195360195)
    r35 = convert(T, -7.214545454545455)
    r34 = convert(T, 20.676923076923078)
    r33 = convert(T, -20.31142191142191)
    r32 = convert(T, 7.14965034965035)
    r45 = convert(T, 7.866666666666666)
    r44 = convert(T, -18.78205128205128)
    r43 = convert(T, 13.508547008547009)
    r42 = convert(T, -2.3653846153846154)
    r55 = convert(T, 2)
    r54 = convert(T, -5.294871794871795)
    r53 = convert(T, 4.534188034188034)
    r52 = convert(T, -1.2115384615384615)
    r65 = convert(T, -1.4924630924630924)
    r64 = convert(T, 1.5785667324128863)
    r63 = convert(T, 1.1958838881915805)
    r62 = convert(T, -1.219801565955412)
    r75 = convert(T, -7.051721611721612)
    r74 = convert(T, 16.273203719357564)
    r73 = convert(T, -11.978886071193763)
    r72 = convert(T, 3.0512256973795435)
    r85 = convert(T, 4)
    r84 = convert(T, -8.384615384615385)
    r83 = convert(T, 5.769230769230769)
    r82 = convert(T, -1.3846153846153846)
    OwrenZen5ConstantCache(a21, a31, a32, a41, a42, a51, a52, a53,
                           a54, a61, a62, a63, a64, a65, a71, a72, a73,
                           a74, a75, a76, a81, a83, a84, a85, a86, a87,
                           c1, c2, c3, c4, c5, c6, btilde1, btilde3, btilde4, btilde5,
                           btilde6, btilde7,
                           r15, r14, r13, r12, r35, r34, r33, r32, r45, r44,
                           r43, r42, r55, r54, r53, r52, r65, r64, r63,
                           r62, r75, r74, r73, r72, r85, r84, r83, r82)
end

function OwrenZen5ConstantCache(T, T2)
    a21 = convert(T, 1 // 6)
    a31 = convert(T, 1 // 16)
    a32 = convert(T, 3 // 16)
    a41 = convert(T, 1 // 4)
    a42 = convert(T, -3 // 4)
    a43 = convert(T, 1)
    a51 = convert(T, -3 // 4)
    a52 = convert(T, 15 // 4)
    a53 = convert(T, -3)
    a54 = convert(T, 1 // 2)
    a61 = convert(T, 369 // 1372)
    a62 = convert(T, -243 // 343)
    a63 = convert(T, 297 // 343)
    a64 = convert(T, 1485 // 9604)
    a65 = convert(T, 297 // 4802)
    a71 = convert(T, -133 // 4512)
    a72 = convert(T, 1113 // 6016)
    a73 = convert(T, 7945 // 16544)
    a74 = convert(T, -12845 // 24064)
    a75 = convert(T, -315 // 24064)
    a76 = convert(T, 156065 // 198528)
    a81 = convert(T, 83 // 945)
    a83 = convert(T, 248 // 825)
    a84 = convert(T, 41 // 180)
    a85 = convert(T, 1 // 36)
    a86 = convert(T, 2401 // 38610)
    a87 = convert(T, 6016 // 20475)
    c1 = convert(T2, 1 // 6)
    c2 = convert(T2, 1 // 4)
    c3 = convert(T2, 1 // 2)
    c4 = convert(T2, 1 // 2)
    c5 = convert(T2, 9 // 14)
    c6 = convert(T2, 7 // 8)
    # b1 = convert(T,-1//9)
    # b3 = convert(T,40//33)
    # b4 = convert(T,-7//4)
    # b5 = convert(T,-1//12)
    # b6 = convert(T,343//198)
    btilde1 = convert(T, -188 // 945)
    btilde3 = convert(T, 752 // 825)
    btilde4 = convert(T, -89 // 45)
    btilde5 = convert(T, -1 // 9)
    btilde6 = convert(T, 32242 // 19305)
    btilde7 = convert(T, -6016 // 20475)
    r15 = convert(T, 596 // 315)
    r14 = convert(T, -4969 // 819)
    r13 = convert(T, 17893 // 2457)
    r12 = convert(T, -3292 // 819)
    r35 = convert(T, -1984 // 275)
    r34 = convert(T, 1344 // 65)
    r33 = convert(T, -43568 // 2145)
    r32 = convert(T, 5112 // 715)
    r45 = convert(T, 118 // 15)
    r44 = convert(T, -1465 // 78)
    r43 = convert(T, 3161 // 234)
    r42 = convert(T, -123 // 52)
    r55 = convert(T, 2)
    r54 = convert(T, -413 // 78)
    r53 = convert(T, 1061 // 234)
    r52 = convert(T, -63 // 52)
    r65 = convert(T, -9604 // 6435)
    r64 = convert(T, 2401 // 1521)
    r63 = convert(T, 60025 // 50193)
    r62 = convert(T, -40817 // 33462)
    r75 = convert(T, -48128 // 6825)
    r74 = convert(T, 96256 // 5915)
    r73 = convert(T, -637696 // 53235)
    r72 = convert(T, 18048 // 5915)
    r85 = convert(T, 4)
    r84 = convert(T, -109 // 13)
    r83 = convert(T, 75 // 13)
    r82 = convert(T, -18 // 13)
    OwrenZen5ConstantCache(a21, a31, a32, a41, a42, a51, a52, a53,
                           a54, a61, a62, a63, a64, a65, a71, a72, a73,
                           a74, a75, a76, a81, a83, a84, a85, a86, a87,
                           c1, c2, c3, c4, c5, c6, btilde1, btilde3, btilde4, btilde5,
                           btilde6, btilde7,
                           r15, r14, r13, r12, r35, r34, r33, r32, r45, r44,
                           r43, r42, r55, r54, r53, r52, r65, r64, r63,
                           r62, r75, r74, r73, r72, r85, r84, r83, r82)
end

struct Tsit5ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
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
    a72::T
    a73::T
    a74::T
    a75::T
    a76::T
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    btilde5::T
    btilde6::T
    btilde7::T
    r11::T
    r12::T
    r13::T
    r14::T
    r22::T
    r23::T
    r24::T
    r32::T
    r33::T
    r34::T
    r42::T
    r43::T
    r44::T
    r52::T
    r53::T
    r54::T
    r62::T
    r63::T
    r64::T
    r72::T
    r73::T
    r74::T
end

function Tsit5ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    c1 = convert(T2, 0.161)
    c2 = convert(T2, 0.327)
    c3 = convert(T2, 0.9)
    c4 = convert(T2, 0.9800255409045097)
    c5 = convert(T2, 1)
    c6 = convert(T2, 1)
    a21 = convert(T, 0.161)
    a31 = convert(T, -0.008480655492356989)
    a32 = convert(T, 0.335480655492357)
    a41 = convert(T, 2.8971530571054935)
    a42 = convert(T, -6.359448489975075)
    a43 = convert(T, 4.3622954328695815)
    a51 = convert(T, 5.325864828439257)
    a52 = convert(T, -11.748883564062828)
    a53 = convert(T, 7.4955393428898365)
    a54 = convert(T, -0.09249506636175525)
    a61 = convert(T, 5.86145544294642)
    a62 = convert(T, -12.92096931784711)
    a63 = convert(T, 8.159367898576159)
    a64 = convert(T, -0.071584973281401)
    a65 = convert(T, -0.028269050394068383)
    a71 = convert(T, 0.09646076681806523)
    a72 = convert(T, 0.01)
    a73 = convert(T, 0.4798896504144996)
    a74 = convert(T, 1.379008574103742)
    a75 = convert(T, -3.290069515436081)
    a76 = convert(T, 2.324710524099774)
    # b1 =       convert(T,0.09468075576583945)
    # b2 =       convert(T,0.009183565540343254)
    # b3 =       convert(T,0.4877705284247616)
    # b4 =       convert(T,1.234297566930479)
    # b5 =       convert(T,-2.7077123499835256)
    # b6 =       convert(T,1.866628418170587)
    # b7 =       convert(T,0.015151515151515152)
    btilde1 = convert(T, -0.00178001105222577714)
    btilde2 = convert(T, -0.0008164344596567469)
    btilde3 = convert(T, 0.007880878010261995)
    btilde4 = convert(T, -0.1447110071732629)
    btilde5 = convert(T, 0.5823571654525552)
    btilde6 = convert(T, -0.45808210592918697)
    btilde7 = convert(T, 0.015151515151515152)

    r11, r12, r13, r14, r22, r23, r24, r32, r33, r34, r42, r43, r44, r52, r53, r54, r62, r63, r64, r72, r73, r74 = Tsit5Interp(T)
    Tsit5ConstantCache(c1, c2, c3, c4, c5, c6, a21, a31, a32, a41, a42, a43, a51, a52, a53,
                       a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, btilde1,
                       btilde2, btilde3, btilde4, btilde5, btilde6, btilde7, r11, r12, r13,
                       r14, r22, r23, r24, r32, r33, r34, r42, r43, r44, r52, r53, r54, r62,
                       r63, r64, r72, r73, r74)
end

function Tsit5ConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 161 // 1000)
    c2 = convert(T2, 327 // 1000)
    c3 = convert(T2, 9 // 10)
    c4 = convert(T2,
                 big".9800255409045096857298102862870245954942137979563024768854764293221195950761080302604")
    c5 = convert(T2, 1)
    c6 = convert(T2, 1)
    a21 = convert(T, 161 // 1000)
    a31 = convert(T,
                  big"-.8480655492356988544426874250230774675121177393430391537369234245294192976164141156943e-2")
    a32 = convert(T,
                  big".3354806554923569885444268742502307746751211773934303915373692342452941929761641411569")
    a41 = convert(T,
                  big"2.897153057105493432130432594192938764924887287701866490314866693455023795137503079289")
    a42 = convert(T,
                  big"-6.359448489975074843148159912383825625952700647415626703305928850207288721235210244366")
    a43 = convert(T,
                  big"4.362295432869581411017727318190886861027813359713760212991062156752264926097707165077")
    a51 = convert(T,
                  big"5.325864828439256604428877920840511317836476253097040101202360397727981648835607691791")
    a52 = convert(T,
                  big"-11.74888356406282787774717033978577296188744178259862899288666928009020615663593781589")
    a53 = convert(T,
                  big"7.495539342889836208304604784564358155658679161518186721010132816213648793440552049753")
    a54 = convert(T,
                  big"-.9249506636175524925650207933207191611349983406029535244034750452930469056411389539635e-1")
    a61 = convert(T,
                  big"5.861455442946420028659251486982647890394337666164814434818157239052507339770711679748")
    a62 = convert(T,
                  big"-12.92096931784710929170611868178335939541780751955743459166312250439928519268343184452")
    a63 = convert(T,
                  big"8.159367898576158643180400794539253485181918321135053305748355423955009222648673734986")
    a64 = convert(T,
                  big"-.7158497328140099722453054252582973869127213147363544882721139659546372402303777878835e-1")
    a65 = convert(T,
                  big"-.2826905039406838290900305721271224146717633626879770007617876201276764571291579142206e-1")
    a71 = convert(T,
                  big".9646076681806522951816731316512876333711995238157997181903319145764851595234062815396e-1")
    a72 = convert(T, 1 // 100)
    a73 = convert(T,
                  big".4798896504144995747752495322905965199130404621990332488332634944254542060153074523509")
    a74 = convert(T,
                  big"1.379008574103741893192274821856872770756462643091360525934940067397245698027561293331")
    a75 = convert(T,
                  big"-3.290069515436080679901047585711363850115683290894936158531296799594813811049925401677")
    a76 = convert(T,
                  big"2.324710524099773982415355918398765796109060233222962411944060046314465391054716027841")
    # b1 =       convert(T,big".9468075576583945807478876255758922856117527357724631226139574065785592789071067303271e-1")
    # b2 =       convert(T,big".9183565540343253096776363936645313759813746240984095238905939532922955247253608687270e-2")
    # b3 =       convert(T,big".4877705284247615707855642599631228241516691959761363774365216240304071651579571959813")
    # b4 =       convert(T,big"1.234297566930478985655109673884237654035539930748192848315425833500484878378061439761")
    # b5 =       convert(T,big"-2.707712349983525454881109975059321670689605166938197378763992255714444407154902012702")
    # b6 =       convert(T,big"1.866628418170587035753719399566211498666255505244122593996591602841258328965767580089")
    # b7 =       convert(T,1//66)
    btilde1 = convert(T,
                      big"-1.780011052225771443378550607539534775944678804333659557637450799792588061629796e-03")
    btilde2 = convert(T,
                      big"-8.164344596567469032236360633546862401862537590159047610940604670770447527463931e-04")
    btilde3 = convert(T,
                      big"7.880878010261996010314727672526304238628733777103128603258129604952959142646516e-03")
    btilde4 = convert(T,
                      big"-1.44711007173262907537165147972635116720922712343167677619514233896760819649515e-01")
    btilde5 = convert(T,
                      big"5.823571654525552250199376106520421794260781239567387797673045438803694038950012e-01")
    btilde6 = convert(T,
                      big"-4.580821059291869466616365188325542974428047279788398179474684434732070620889539e-01")
    btilde7 = convert(T, 1 // 66)

    r11, r12, r13, r14, r22, r23, r24, r32, r33, r34, r42, r43, r44, r52, r53, r54, r62, r63, r64, r72, r73, r74 = Tsit5Interp(T)
    Tsit5ConstantCache(c1, c2, c3, c4, c5, c6, a21, a31, a32, a41, a42, a43, a51, a52, a53,
                       a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, btilde1,
                       btilde2, btilde3, btilde4, btilde5, btilde6, btilde7, r11, r12, r13,
                       r14, r22, r23, r24, r32, r33, r34, r42, r43, r44, r52, r53, r54, r62,
                       r63, r64, r72, r73, r74)
end

"""
Coefficients for the polynomial
bᵢΘ = ri1*Θ + ri2*Θ^2 + ri3*Θ^3 + ...

These are the coefficients of the expanded form of the polynomials from

Runge–Kutta pairs of order 5(4) satisfying only the first column
simplifying assumption

Ch. Tsitouras
"""
function Tsit5Interp(T::Type{<:CompiledFloats})
    r11 = convert(T, 1.0)
    r12 = convert(T, -2.763706197274826)
    r13 = convert(T, 2.9132554618219126)
    r14 = convert(T, -1.0530884977290216)

    r22 = convert(T, 0.13169999999999998)
    r23 = convert(T, -0.2234)
    r24 = convert(T, 0.1017)

    r32 = convert(T, 3.9302962368947516)
    r33 = convert(T, -5.941033872131505)
    r34 = convert(T, 2.490627285651253)

    r42 = convert(T, -12.411077166933676)
    r43 = convert(T, 30.33818863028232)
    r44 = convert(T, -16.548102889244902)

    r52 = convert(T, 37.50931341651104)
    r53 = convert(T, -88.1789048947664)
    r54 = convert(T, 47.37952196281928)

    r62 = convert(T, -27.896526289197286)
    r63 = convert(T, 65.09189467479366)
    r64 = convert(T, -34.87065786149661)

    r72 = convert(T, 1.5)
    r73 = convert(T, -4)
    r74 = convert(T, 2.5)

    return r11, r12, r13, r14, r22, r23, r24, r32, r33, r34, r42, r43, r44, r52, r53, r54,
           r62, r63, r64, r72, r73, r74
end

"""
Coefficients for the polynomial
bᵢΘ = ri1*Θ + ri2*Θ^2 + ri3*Θ^3 + ...

These are the coefficients of the expanded form of the polynomials from

Runge–Kutta pairs of order 5(4) satisfying only the first column
simplifying assumption

Ch. Tsitouras
"""
function Tsit5Interp(T::Type)
    r11 = convert(T, big"0.999999999999999974283372471559910888475488471328")
    r12 = convert(T, big"-2.763706197274825911336735930481400260916070804192")
    r13 = convert(T, big"2.91325546182191274375068099306808")
    r14 = convert(T, -1.0530884977290216)

    r22 = convert(T, big"0.13169999999999999727")
    r23 = convert(T, big"-0.22339999999999999818")
    r24 = convert(T, 0.1017)

    r32 = convert(T, big"3.93029623689475152850687446709813398")
    r33 = convert(T, big"-5.94103387213150473470249202589458001")
    r34 = convert(T, big"2.490627285651252793")

    r42 = convert(T, big"-12.411077166933676983734381540685453484102414134010752")
    r43 = convert(T, big"30.3381886302823215981729903691836576")
    r44 = convert(T, big"-16.54810288924490272")

    r52 = convert(T, big"37.50931341651103919496903965334519631242339792120440212")
    r53 = convert(T, big"-88.1789048947664011014276693541209817")
    r54 = convert(T, big"47.37952196281928122")

    r62 = convert(T, big"-27.896526289197287805948263144598643896")
    r63 = convert(T, big"65.09189467479367152629021928716553658")
    r64 = convert(T, big"-34.87065786149660974")

    r72 = convert(T, 1.5)
    r73 = convert(T, -4)
    r74 = convert(T, 2.5)

    return r11, r12, r13, r14, r22, r23, r24, r32, r33, r34, r42, r43, r44, r52, r53, r54,
           r62, r63, r64, r72, r73, r74
end

struct BS5ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
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
    a71::T
    a72::T
    a73::T
    a74::T
    a75::T
    a76::T
    a81::T
    a83::T
    a84::T
    a85::T
    a86::T
    a87::T
    bhat1::T
    bhat3::T
    bhat4::T
    bhat5::T
    bhat6::T
    btilde1::T
    btilde3::T
    btilde4::T
    btilde5::T
    btilde6::T
    btilde7::T
    btilde8::T
    c6::T2
    c7::T2
    c8::T2
    a91::T
    a92::T
    a93::T
    a94::T
    a95::T
    a96::T
    a97::T
    a98::T
    a101::T
    a102::T
    a103::T
    a104::T
    a105::T
    a106::T
    a107::T
    a108::T
    a109::T
    a111::T
    a112::T
    a113::T
    a114::T
    a115::T
    a116::T
    a117::T
    a118::T
    a119::T
    a1110::T
    r016::T
    r015::T
    r014::T
    r013::T
    r012::T
    r036::T
    r035::T
    r034::T
    r033::T
    r032::T
    r046::T
    r045::T
    r044::T
    r043::T
    r042::T
    r056::T
    r055::T
    r054::T
    r053::T
    r052::T
    r066::T
    r065::T
    r064::T
    r063::T
    r062::T
    r076::T
    r075::T
    r074::T
    r073::T
    r072::T
    r086::T
    r085::T
    r084::T
    r083::T
    r082::T
    r096::T
    r095::T
    r094::T
    r093::T
    r106::T
    r105::T
    r104::T
    r103::T
    r102::T
    r116::T
    r115::T
    r114::T
    r113::T
    r112::T
end

"""

An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28
"""
function BS5ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    c1 = convert(T2, 0.16666666666666666)
    c2 = convert(T2, 0.2222222222222222)
    c3 = convert(T2, 0.42857142857142855)
    c4 = convert(T2, 0.6666666666666666)
    c5 = convert(T2, 0.75)

    a21 = convert(T, 0.16666666666666666)
    a31 = convert(T, 0.07407407407407407)
    a32 = convert(T, 0.14814814814814814)
    a41 = convert(T, 0.13338192419825073)
    a42 = convert(T, -0.47230320699708456)
    a43 = convert(T, 0.7674927113702624)
    a51 = convert(T, 0.22895622895622897)
    a52 = convert(T, -0.36363636363636365)
    a53 = convert(T, 0.2937062937062937)
    a54 = convert(T, 0.5076405076405076)
    a61 = convert(T, 0.026500355113636364)
    a62 = convert(T, 0.23011363636363635)
    a63 = convert(T, 0.10772747760052448)
    a64 = convert(T, 0.1602190777972028)
    a65 = convert(T, 0.225439453125)
    a71 = convert(T, 0.18159821692916506)
    a72 = convert(T, -0.38707982536247293)
    a73 = convert(T, 0.41288268560074054)
    a74 = convert(T, 0.6409913191253497)
    a75 = convert(T, -1.0121439383514517)
    a76 = convert(T, 1.1637515420586695)
    a81 = convert(T, 0.07279265873015874)
    a83 = convert(T, 0.28662437773692473)
    a84 = convert(T, 0.19513621794871794)
    a85 = convert(T, 0.008638392857142857)
    a86 = convert(T, 0.3595655806182122)
    a87 = convert(T, 0.07724277210884353)

    bhat1 = convert(T, -0.00234375)
    bhat3 = convert(T, 0.0103760754048583)
    bhat4 = convert(T, -0.016490384615384615)
    bhat5 = convert(T, 0.018984375)
    bhat6 = convert(T, -0.010526315789473684)
    #  btilde1=convert(T,0.07084476451760402)
    #  btilde2=convert(T,0)
    #  btilde3=convert(T,0.2956730769230769)
    #  btilde4=convert(T,0.17965747482208388)
    #  btilde5=convert(T,0.029861111111111113)
    #  btilde6=convert(T,0.3462886755067825)
    #  btilde7=convert(T,0.07176240133870539)
    btilde1 = convert(T, -0.0019478942125547064)
    btilde3 = convert(T, 0.009048699186152193)
    btilde4 = convert(T, -0.015478743126634073)
    btilde5 = convert(T, 0.021222718253968254)
    btilde6 = convert(T, -0.013276905111429694)
    btilde7 = convert(T, -0.005480370770138146)
    btilde8 = convert(T, 0.005912495780636172)
    c6, c7, c8, a91, a92, a93, a94, a95, a96, a97, a98, a101, a102, a103, a104, a105, a106, a107, a108, a109, a111, a112, a113, a114, a115, a116, a117, a118, a119, a1110 = BS5Interp(T,
                                                                                                                                                                                      T2)
    r016, r015, r014, r013, r012, r036, r035, r034, r033, r032, r046, r045, r044, r043, r042, r056, r055, r054, r053, r052, r066, r065, r064, r063, r062, r076, r075, r074, r073, r072, r086, r085, r084, r083, r082, r096, r095, r094, r093, r106, r105, r104, r103, r102, r116, r115, r114, r113, r112 = BS5Interp_polyweights(T)
    BS5ConstantCache(c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54,
                     a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a83, a84,
                     a85, a86, a87, bhat1, bhat3, bhat4, bhat5, bhat6, btilde1, btilde3,
                     btilde4, btilde5, btilde6, btilde7, btilde8, c6, c7, c8, a91, a92, a93,
                     a94, a95, a96, a97, a98, a101, a102, a103, a104, a105, a106, a107,
                     a108, a109, a111, a112, a113, a114, a115, a116, a117, a118, a119,
                     a1110, r016, r015, r014, r013, r012, r036, r035, r034, r033, r032,
                     r046, r045, r044, r043, r042, r056, r055, r054, r053, r052, r066, r065,
                     r064, r063, r062, r076, r075, r074, r073, r072, r086, r085, r084, r083,
                     r082, r096, r095, r094, r093, r106, r105, r104, r103, r102, r116, r115,
                     r114, r113, r112)
end

"""

An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28
"""
function BS5ConstantCache(T::Type, T2::Type)
    c1 = convert(T2, 1 // 6)
    c2 = convert(T2, 2 // 9)
    c3 = convert(T2, 3 // 7)
    c4 = convert(T2, 2 // 3)
    c5 = convert(T2, 3 // 4)

    a21 = convert(T, 1 // 6)
    a31 = convert(T, 2 // 27)
    a32 = convert(T, 4 // 27)
    a41 = convert(T, 183 // 1372)
    a42 = convert(T, -162 // 343)
    a43 = convert(T, 1053 // 1372)
    a51 = convert(T, 68 // 297)
    a52 = convert(T, -4 // 11)
    a53 = convert(T, 42 // 143)
    a54 = convert(T, 1960 // 3861)
    a61 = convert(T, 597 // 22528)
    a62 = convert(T, 81 // 352)
    a63 = convert(T, 63099 // 585728)
    a64 = convert(T, 58653 // 366080)
    a65 = convert(T, 4617 // 20480)
    a71 = convert(T, 174197 // 959244)
    a72 = convert(T, -30942 // 79937)
    a73 = convert(T, 8152137 // 19744439)
    a74 = convert(T, 666106 // 1039181)
    a75 = convert(T, -29421 // 29068)
    a76 = convert(T, 482048 // 414219)
    a81 = convert(T, 587 // 8064)
    a83 = convert(T, 4440339 // 15491840)
    a84 = convert(T, 24353 // 124800)
    a85 = convert(T, 387 // 44800)
    a86 = convert(T, 2152 // 5985)
    a87 = convert(T, 7267 // 94080)

    bhat1 = convert(T, -3 // 1280)
    bhat3 = convert(T, 6561 // 632320)
    bhat4 = convert(T, -343 // 20800)
    bhat5 = convert(T, 243 // 12800)
    bhat6 = convert(T, -1 // 95)
    #  btilde1=convert(T,2479//34992)
    #  btilde2=convert(T,0)
    #  btilde3=convert(T,123//416)
    #  btilde4=convert(T,612941//3411720)
    #  btilde5=convert(T,43//1440)
    #  btilde6=convert(T,2272//6561)
    #  btilde7=convert(T,79937//1113912)
    btilde1 = convert(T, -3817 // 1959552)
    btilde3 = convert(T, 140181 // 15491840)
    btilde4 = convert(T, -4224731 // 272937600)
    btilde5 = convert(T, 8557 // 403200)
    btilde6 = convert(T, -57928 // 4363065)
    btilde7 = convert(T, -23930231 // 4366535040)
    btilde8 = convert(T, 3293 // 556956)
    c6, c7, c8, a91, a92, a93, a94, a95, a96, a97, a98, a101, a102, a103, a104, a105, a106, a107, a108, a109, a111, a112, a113, a114, a115, a116, a117, a118, a119, a1110 = BS5Interp(T,
                                                                                                                                                                                      T2)
    r016, r015, r014, r013, r012, r036, r035, r034, r033, r032, r046, r045, r044, r043, r042, r056, r055, r054, r053, r052, r066, r065, r064, r063, r062, r076, r075, r074, r073, r072, r086, r085, r084, r083, r082, r096, r095, r094, r093, r106, r105, r104, r103, r102, r116, r115, r114, r113, r112 = BS5Interp_polyweights(T)
    BS5ConstantCache(c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54,
                     a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a83, a84,
                     a85, a86, a87, bhat1, bhat3, bhat4, bhat5, bhat6, btilde1, btilde3,
                     btilde4, btilde5, btilde6, btilde7, btilde8, c6, c7, c8, a91, a92, a93,
                     a94, a95, a96, a97, a98, a101, a102, a103, a104, a105, a106, a107,
                     a108, a109, a111, a112, a113, a114, a115, a116, a117, a118, a119,
                     a1110, r016, r015, r014, r013, r012, r036, r035, r034, r033, r032,
                     r046, r045, r044, r043, r042, r056, r055, r054, r053, r052, r066, r065,
                     r064, r063, r062, r076, r075, r074, r073, r072, r086, r085, r084, r083,
                     r082, r096, r095, r094, r093, r106, r105, r104, r103, r102, r116, r115,
                     r114, r113, r112)
end

"""

An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28


Used in the lazy construction of the dense output

k9, k10, k11 are not computed until called in the dense routine
"""
function BS5Interp(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    c6 = convert(T2, 0.5)
    c7 = convert(T2, 0.8333333333333334)
    c8 = convert(T2, 0.1111111111111111)
    a91 = convert(T, 0.07405598958333333)
    a92 = convert(T, 0)
    a93 = convert(T, 0.28964485093442743)
    a94 = convert(T, 0.12839214966168092)
    a95 = convert(T, -0.003779296875)
    a96 = convert(T, 0.014230019493177388)
    a97 = convert(T, -0.03379371279761905)
    a98 = convert(T, 0.03125)
    a101 = convert(T, -0.06358724036162344)
    a102 = convert(T, 0.5742461924818869)
    a103 = convert(T, -0.06365063007249953)
    a104 = convert(T, 0.043159777438314964)
    a105 = convert(T, 0.8370112883898733)
    a106 = convert(T, -0.34045447246719235)
    a107 = convert(T, 0.04926503818334922)
    a108 = convert(T, -0.006882677669165967)
    a109 = convert(T, -0.19577394258960973)
    a111 = convert(T, 0.0636090772400987)
    a112 = convert(T, 0.01057854182854183)
    a113 = convert(T, 0.06600100945670531)
    a114 = convert(T, 0.02048391555358402)
    a115 = convert(T, 0.003682270330219549)
    a116 = convert(T, 0.155258632271002)
    a117 = convert(T, -0.08509702513818027)
    a118 = convert(T, 0.1)
    a119 = convert(T, -0.1)
    a1110 = convert(T, -0.12340531043086005)

    return c6, c7, c8, a91, a92, a93, a94, a95, a96, a97, a98, a101, a102, a103, a104, a105,
           a106, a107, a108, a109, a111, a112, a113, a114, a115, a116, a117, a118, a119,
           a1110
end

"""

An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28


Used in the lazy construction of the dense output

k9, k10, k11 are not computed until called in the dense routine
"""
function BS5Interp(T::Type, T2::Type)
    c6 = convert(T2, 1 // 2)
    c7 = convert(T2, 5 // 6)
    c8 = convert(T2, 1 // 9)
    a91 = convert(T, 455 // 6144)
    a92 = convert(T, 0)
    a93 = convert(T, 10256301 // 35409920)
    a94 = convert(T, 2307361 // 17971200)
    a95 = convert(T, -387 // 102400)
    a96 = convert(T, 73 // 5130)
    a97 = convert(T, -7267 // 215040)
    a98 = convert(T, 1 // 32)
    a101 = convert(T, -837888343715 // 13176988637184)
    a102 = convert(T, 30409415 // 52955362)
    a103 = convert(T, -48321525963 // 759168069632)
    a104 = convert(T, 8530738453321 // 197654829557760)
    a105 = convert(T, 1361640523001 // 1626788720640)
    a106 = convert(T, -13143060689 // 38604458898)
    a107 = convert(T, 18700221969 // 379584034816)
    a108 = convert(T, -5831595 // 847285792)
    a109 = convert(T, -5183640 // 26477681)
    a111 = convert(T, 98719073263 // 1551965184000)
    a112 = convert(T, 1307 // 123552)
    a113 = convert(T, 4632066559387 // 70181753241600)
    a114 = convert(T, 7828594302389 // 382182512025600)
    a115 = convert(T, 40763687 // 11070259200)
    a116 = convert(T, 34872732407 // 224610586200)
    a117 = convert(T, -2561897 // 30105600)
    a118 = convert(T, 1 // 10)
    a119 = convert(T, -1 // 10)
    a1110 = convert(T, -1403317093 // 11371610250)

    return c6, c7, c8, a91, a92, a93, a94, a95, a96, a97, a98, a101, a102, a103, a104, a105,
           a106, a107, a108, a109, a111, a112, a113, a114, a115, a116, a117, a118, a119,
           a1110
end

"""
Coefficients for the polynomial
bᵢΘ = ri1*Θ + ri2*Θ^2 + ri3*Θ^3 + ...

These coefficients are taken from RKSuite

Note that RKSuite has an error: r081 should be 0
and r011 should be 1. This is pretty easy to spot
since the first order interpolation is linear from y₀.
"""
function BS5Interp_polyweights(T::Type{<:CompiledFloats})
    r016 = convert(T, -11.547607240534195)
    r015 = convert(T, 36.89579791207878)
    r014 = convert(T, -45.470343197183475)
    r013 = convert(T, 27.298136302807084)
    r012 = convert(T, -8.103191118438032)

    r036 = convert(T, -27.245925284262768)
    r035 = convert(T, 74.85879078710211)
    r034 = convert(T, -72.2759709558427)
    r033 = convert(T, 28.38602193195626)
    r032 = convert(T, -3.4362921012159933)

    r046 = convert(T, -14.307769126529058)
    r045 = convert(T, 38.240038148817945)
    r044 = convert(T, -35.469854845413806)
    r043 = convert(T, 13.060399314592578)
    r042 = convert(T, -1.3276772735189397)

    r056 = convert(T, -0.7296779223862557)
    r055 = convert(T, 1.9817123385873385)
    r054 = convert(T, -1.964240395021645)
    r053 = convert(T, 0.8847786781120115)
    r052 = convert(T, -0.16393430643430643)

    r066 = convert(T, 19.121088881250603)
    r065 = convert(T, -65.9928405785889)
    r064 = convert(T, 77.4690381021765)
    r063 = convert(T, -34.163041154825144)
    r062 = convert(T, 3.9253203306051474)

    r076 = convert(T, -14.676126700680273)
    r075 = convert(T, 42.17455357142857)
    r074 = convert(T, -42.4062818877551)
    r073 = convert(T, 16.83892431972789)
    r072 = convert(T, -1.853826530612245)

    r086 = convert(T, 19.453125)
    r085 = convert(T, -54.359375)
    r084 = convert(T, 53.671875)
    r083 = convert(T, -21.078125)
    r082 = convert(T, 2.3125)

    r096 = convert(T, 18.333333333333332)
    r095 = convert(T, -39)
    r094 = convert(T, 23)
    r093 = convert(T, -2.3333333333333335)

    r106 = convert(T, -23.400440940191388)
    r105 = convert(T, 70.20132282057416)
    r104 = convert(T, -73.55422182095978)
    r103 = convert(T, 30.106238940962648)
    r102 = convert(T, -3.3528990003856314)

    r116 = convert(T, 35)
    r115 = convert(T, -105)
    r114 = convert(T, 117)
    r113 = convert(T, -59)
    r112 = convert(T, 12)

    return r016, r015, r014, r013, r012, r036, r035, r034, r033, r032, r046, r045, r044,
           r043, r042, r056, r055, r054, r053, r052, r066, r065, r064, r063, r062, r076,
           r075, r074, r073, r072, r086, r085, r084, r083, r082, r096, r095, r094, r093,
           r106, r105, r104, r103, r102, r116, r115, r114, r113, r112
end

"""
Coefficients for the polynomial
bᵢΘ = ri1*Θ + ri2*Θ^2 + ri3*Θ^3 + ...

These coefficients are taken from RKSuite

Note that RKSuite has an error: r081 should be 0
and r011 should be 1. This is pretty easy to spot
since the first order interpolation is linear from y₀.
"""
function BS5Interp_polyweights(T::Type)
    r016 = convert(T, -12134338393 // 1050809760)
    r015 = convert(T, 12923488183 // 350269920)
    r014 = convert(T, -2722545893 // 59875200)
    r013 = convert(T, 35856435071 // 1313512200)
    r012 = convert(T, -3547880131 // 437837400)

    r036 = convert(T, -33197340367 // 1218433216)
    r035 = convert(T, 65150312289 // 870309440)
    r034 = convert(T, -27096444225 // 374902528)
    r033 = convert(T, 4323308999 // 152304152)
    r032 = convert(T, -1046723109 // 304608304)

    r046 = convert(T, -284800997201 // 19905339168)
    r045 = convert(T, 6343174409579 // 165877826400)
    r044 = convert(T, -201150852119 // 5671036800)
    r043 = convert(T, 3249645975331 // 248816739600)
    r042 = convert(T, -55058055073 // 41469456600)

    r056 = convert(T, -540919 // 741312)
    r055 = convert(T, 85695583 // 43243200)
    r054 = convert(T, -2903933 // 1478400)
    r053 = convert(T, 3586937 // 4054050)
    r052 = convert(T, -1772261 // 10810800)

    r066 = convert(T, 7157998304 // 374350977)
    r065 = convert(T, -41174140576 // 623918295)
    r064 = convert(T, 413114104 // 5332635)
    r063 = convert(T, -9134977024 // 267393555)
    r062 = convert(T, 2449079168 // 623918295)

    r076 = convert(T, -138073 // 9408)
    r075 = convert(T, 94471 // 2240)
    r074 = convert(T, -1329861 // 31360)
    r073 = convert(T, 792103 // 47040)
    r072 = convert(T, -7267 // 3920)

    r086 = convert(T, 1245 // 64)
    r085 = convert(T, -3479 // 64)
    r084 = convert(T, 3435 // 64)
    r083 = convert(T, -1349 // 64)
    r082 = convert(T, 37 // 16)

    r096 = convert(T, 55 // 3)
    r095 = convert(T, -39)
    r094 = convert(T, 23)
    r093 = convert(T, -7 // 3)

    r106 = convert(T, -1774004627 // 75810735)
    r105 = convert(T, 1774004627 // 25270245)
    r104 = convert(T, -26477681 // 359975)
    r103 = convert(T, 11411880511 // 379053675)
    r102 = convert(T, -423642896 // 126351225)

    r116 = convert(T, 35)
    r115 = convert(T, -105)
    r114 = convert(T, 117)
    r113 = convert(T, -59)
    r112 = convert(T, 12)

    return r016, r015, r014, r013, r012, r036, r035, r034, r033, r032, r046, r045, r044,
           r043, r042, r056, r055, r054, r053, r052, r066, r065, r064, r063, r062, r076,
           r075, r074, r073, r072, r086, r085, r084, r083, r082, r096, r095, r094, r093,
           r106, r105, r104, r103, r102, r116, r115, r114, r113, r112
end

struct DP5ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
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
    btilde1::T
    btilde3::T
    btilde4::T
    btilde5::T
    btilde6::T
    btilde7::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2
    d1::T
    d3::T
    d4::T
    d5::T
    d6::T
    d7::T
end

function DP5ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    a21 = convert(T, 0.2)
    a31 = convert(T, 0.075)
    a32 = convert(T, 0.225)
    a41 = convert(T, 0.9777777777777777)
    a42 = convert(T, -3.7333333333333334)
    a43 = convert(T, 3.5555555555555554)
    a51 = convert(T, 2.9525986892242035)
    a52 = convert(T, -11.595793324188385)
    a53 = convert(T, 9.822892851699436)
    a54 = convert(T, -0.2908093278463649)
    a61 = convert(T, 2.8462752525252526)
    a62 = convert(T, -10.757575757575758)
    a63 = convert(T, 8.906422717743473)
    a64 = convert(T, 0.2784090909090909)
    a65 = convert(T, -0.2735313036020583)
    a71 = convert(T, 0.09114583333333333)
    a73 = convert(T, 0.44923629829290207)
    a74 = convert(T, 0.6510416666666666)
    a75 = convert(T, -0.322376179245283)
    a76 = convert(T, 0.13095238095238096)
    # b1  = convert(T,0.08991319444444444)
    # b3  = convert(T,0.4534890685834082)
    # b4  = convert(T,0.6140625)
    # b5  = convert(T,-0.2715123820754717)
    # b6  = convert(T,0.08904761904761904)
    # b7  = convert(T,0.025)
    btilde1 = convert(T, -0.0012326388888888888)
    btilde3 = convert(T, 0.0042527702905061394)
    btilde4 = convert(T, -0.03697916666666667)
    btilde5 = convert(T, 0.05086379716981132)
    btilde6 = convert(T, -0.0419047619047619)
    btilde7 = convert(T, 0.025)
    c1 = convert(T2, 0.2)
    c2 = convert(T2, 0.3)
    c3 = convert(T2, 0.8)
    c4 = convert(T2, 0.8888888888888888)
    c5 = convert(T2, 1)
    c6 = convert(T2, 1)
    d1, d3, d4, d5, d6, d7 = DP5_dense_ds(T)
    DP5ConstantCache(a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64,
                     a65, a71, a73, a74, a75, a76, btilde1, btilde3, btilde4, btilde5,
                     btilde6, btilde7, c1, c2, c3, c4, c5, c6, d1, d3, d4, d5, d6, d7)
end

function DP5_dense_ds(T::Type{<:CompiledFloats})
    d1 = convert(T, -1.1270175653862835)
    d3 = convert(T, 2.675424484351598)
    d4 = convert(T, -5.685526961588504)
    d5 = convert(T, 3.5219323679207912)
    d6 = convert(T, -1.7672812570757455)
    d7 = convert(T, 2.382468931778144)
    return d1, d3, d4, d5, d6, d7
end

function DP5ConstantCache(T::Type, T2::Type)
    a21 = convert(T, 1 // 5)
    a31 = convert(T, 3 // 40)
    a32 = convert(T, 9 // 40)
    a41 = convert(T, 44 / 45)
    a42 = convert(T, -56 // 15)
    a43 = convert(T, 32 // 9)
    a51 = convert(T, 19372 // 6561)
    a52 = convert(T, -25360 // 2187)
    a53 = convert(T, 64448 // 6561)
    a54 = convert(T, -212 // 729)
    a61 = convert(T, 9017 // 3168)
    a62 = convert(T, -355 // 33)
    a63 = convert(T, 46732 // 5247)
    a64 = convert(T, 49 // 176)
    a65 = convert(T, -5103 // 18656)
    a71 = convert(T, 35 // 384)
    a73 = convert(T, 500 // 1113)
    a74 = convert(T, 125 // 192)
    a75 = convert(T, -2187 // 6784)
    a76 = convert(T, 11 // 84)
    #  b1  = convert(T,5179//57600)
    #  b3  = convert(T,7571//16695)
    #  b4  = convert(T,393//640)
    #  b5  = convert(T,-92097//339200)
    #  b6  = convert(T,187//2100)
    #  b7  = convert(T,1//40)
    btilde1 = convert(T, -71 // 57600)
    btilde3 = convert(T, 71 // 16695)
    btilde4 = convert(T, -71 // 1920)
    btilde5 = convert(T, 17253 // 339200)
    btilde6 = convert(T, -22 // 525)
    btilde7 = convert(T, 1 // 40)
    c1 = convert(T2, 1 // 5)
    c2 = convert(T2, 3 // 10)
    c3 = convert(T2, 4 // 5)
    c4 = convert(T2, 8 // 9)
    c5 = convert(T2, 1)
    c6 = convert(T2, 1)
    d1, d3, d4, d5, d6, d7 = DP5_dense_ds(T)
    DP5ConstantCache(a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64,
                     a65, a71, a73, a74, a75, a76, btilde1, btilde3, btilde4, btilde5,
                     btilde6, btilde7, c1, c2, c3, c4, c5, c6, d1, d3, d4, d5, d6, d7)
end

function DP5_dense_ds(T)
    d1 = convert(T, -12715105075 // 11282082432)
    d3 = convert(T, 87487479700 // 32700410799)
    d4 = convert(T, -10690763975 // 1880347072)
    d5 = convert(T, 701980252875 // 199316789632)
    d6 = convert(T, -1453857185 // 822651844)
    d7 = convert(T, 69997945 // 29380423)
    return d1, d3, d4, d5, d6, d7
end

#=
function DP5_dense_bs(T)
  b1  = convert(T,5179//57600)
  b3  = convert(T,7571//16695)
  b4  = convert(T,393//640)
  b5  = convert(T,-92097//339200)
  b6  = convert(T,187//2100)
  b7  = convert(T,1//40)
  return b1,b3,b4,b5,b6,b7
end
=#

"""
An Optimized Runge-Kutta method for the solution of Orbital Problems
by Z.A. Anastassi and T.E. Simos
Journal of Computational and Applied Mathematics, Volume 175, Issue 1, 1 March 2005, Pages 1 to 9.
"""
struct Anas5ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
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
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2
    b1::T
    b3::T
    b4::T
    b5::T
    b6::T
end

function Anas5ConstantCache(T::Type{<:CompiledFloats}, T2::Type{<:CompiledFloats})
    a21 = convert(T, 0.1)
    a31 = convert(T, -0.2222222222222222)
    a32 = convert(T, 0.5555555555555556)
    a41 = convert(T, 3.111111111111111)
    a42 = convert(T, -4.444444444444444)
    a43 = convert(T, 2.0)
    a51 = convert(T, -1.409625)
    a52 = convert(T, 2.1375)
    a53 = convert(T, -0.2295)
    a54 = convert(T, 0.401625)
    a61 = convert(T, -4.0)
    a62 = convert(T, 5.0)
    a63 = convert(T, 0.0)
    a64 = convert(T, 0.0)
    a65 = convert(T, 0.0)
    c2 = convert(T2, 0.1)
    c3 = convert(T2, 0.3333333333333333)
    c4 = convert(T2, 0.6666666666666667)
    c5 = convert(T2, 0.9)
    c6 = convert(T2, 1)
    b1 = convert(T, 0.1064814814814815)
    b3 = convert(T, 0.4632352941176471)
    b4 = convert(T, 0.1607142857142857)
    b5 = convert(T, 0.3112356053532524)
    b6 = convert(T, -0.04166666666666667)
    Anas5ConstantCache(a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64,
                       a65, c2, c3, c4, c5, c6, b1, b3, b4, b5, b6)
end

function Anas5ConstantCache(T, T2)
    a21 = convert(T, 1 // 10)
    a31 = convert(T, -2 // 9)
    a32 = convert(T, 5 // 9)
    a41 = convert(T, 28 // 9)
    a42 = convert(T, -40 // 9)
    a43 = convert(T, 2 // 1)
    a51 = convert(T, -11277 // 8000)
    a52 = convert(T, 171 // 80)
    a53 = convert(T, -459 // 2000)
    a54 = convert(T, 3213 // 8000)
    a61 = convert(T, -4 // 1)
    a62 = convert(T, 5 // 1)
    a63 = convert(T, 0 // 1)
    a64 = convert(T, 0 // 1)
    a65 = convert(T, 0 // 1)
    c2 = convert(T2, 1 // 10)
    c3 = convert(T2, 1 // 3)
    c4 = convert(T2, 2 // 3)
    c5 = convert(T2, 9 // 10)
    c6 = convert(T2, 1 // 1)
    b1 = convert(T, 23 // 216)
    b3 = convert(T, 63 // 136)
    b4 = convert(T, 9 // 56)
    b5 = convert(T, 1000 // 3213)
    b6 = convert(T, -1 // 24)
    Anas5ConstantCache(a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64,
                       a65, c2, c3, c4, c5, c6, b1, b3, b4, b5, b6)
end

struct MSRK5ConstantCache{T, T1} <: OrdinaryDiffEqConstantCache
    a21::T
    a31::T
    a32::T
    a41::T
    a43::T
    a51::T
    a53::T
    a54::T
    a61::T
    a63::T
    a64::T
    a65::T
    a71::T
    a73::T
    a74::T
    a75::T
    a76::T
    a81::T
    a83::T
    a84::T
    a85::T
    a86::T
    a87::T

    b1::T # == a_{9i} ∀ i=1:8
    b4::T
    b5::T
    b6::T
    b7::T
    b8::T

    c2::T1
    c3::T1
    c4::T1
    c5::T1
    c6::T1
    c7::T1
    c8::T1
end

# Use rational numbers for testing. Define another function that defines the tab using floats.
function MSRK5ConstantCache(T::Type, T1::Type)
    a21 = T(4 // 45)
    a31 = T(1 // 30)
    a32 = T(1 // 10)
    a41 = T(1 // 20)
    a43 = T(3 // 20)
    a51 = T(1 // 2)
    a53 = T(-15 // 8)
    a54 = T(15 // 8)
    a61 = T(-11 // 135)
    a63 = T(23 // 45)
    a64 = T(-2 // 27)
    a65 = T(8 // 45)
    a71 = T(5 // 108)
    a73 = T(35 // 72)
    a74 = T(-59 // 216)
    a75 = T(-25 // 27)
    a76 = T(3 // 2)
    a81 = T(31 // 128)
    a83 = T(-7563 // 4480)
    a84 = T(233 // 112)
    a85 = T(3461 // 2240)
    a86 = T(-765 // 448)
    a87 = T(153 // 320)
    b1 = T(29 // 456)
    b4 = T(11 // 38)
    b5 = T(2 // 27)
    b6 = T(11 // 40)
    b7 = T(4 // 19)
    b8 = T(224 // 2565)
    c2 = T1(4 // 45)
    c3 = T1(2 // 15)
    c4 = T1(1 // 5)
    c5 = T1(1 // 2)
    c6 = T1(8 // 15)
    c7 = T1(5 // 6)
    c8 = T1(19 // 20)

    MSRK5ConstantCache(a21, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73,
                       a74, a75, a76, a81, a83, a84, a85, a86, a87, b1, b4, b5, b6, b7, b8,
                       c2, c3, c4, c5, c6, c7, c8)
end

struct Stepanov5ConstantCache{T, T1} <: OrdinaryDiffEqConstantCache
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
    # b2::T
    b3::T
    b4::T
    b5::T
    b6::T

    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    btilde5::T
    btilde6::T
    btilde7::T

    c2::T1
    c3::T1
    c4::T1
    c5::T1
    c6::T1
end

# Use rational numbers for testing. Define another function that defines the tab using floats.
function Stepanov5ConstantCache(T::Type, T1::Type)
    a21 = T(1 // 5)
    a31 = T(21 // 338)
    a32 = T(441 // 1690)
    a41 = T(639 // 392)
    a42 = T(-729 // 140)
    a43 = T(1755 // 392)
    a51 = T(4878991 // 1693440)
    a52 = T(-16601 // 1792)
    a53 = T(210067 // 28224)
    a54 = T(-1469 // 17280)
    a61 = T(13759919 // 4230954)
    a62 = T(-2995 // 287)
    a63 = T(507312091 // 61294590)
    a64 = T(-22 // 405)
    a65 = T(-7040 // 180687)
    b1 = T(1441 // 14742)
    # b2 = T(0)
    b3 = T(114244 // 234927)
    b4 = T(118 // 81)
    b5 = T(-12800 // 4407)
    b6 = T(41 // 22)

    btilde1 = T(-1 // 273)
    btilde2 = T(0)
    btilde3 = T(2197 // 174020)
    btilde4 = T(-4 // 15)
    btilde5 = T(1280 // 1469)
    btilde6 = T(-33743 // 52712)
    btilde7 = T(127 // 4792)

    c2 = T1(1 // 5)
    c3 = T1(21 // 65)
    c4 = T1(9 // 10)
    c5 = T1(39 // 40)
    c6 = T1(1 // 1)

    Stepanov5ConstantCache(a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63,
                           a64, a65, b1, b3, b4, b5, b6,
                           btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7,
                           c2, c3, c4, c5, c6)
end
