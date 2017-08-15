struct IRKN3ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  bconst1::T
  bconst2::T
  c1::T2
  a21::T
  b1::T
  b2::T
  bbar1::T
  bbar2::T
end

Base.@pure function IRKN3ConstantCache{T<:CompiledFloats,T2<:CompiledFloats}(::Type{T},::Type{T2})
  bconst1 = T(1.5)
  bconst2 = T(-0.5)
  c1 = T2(0.5)
  a21 = T(0.125)
  b1 = T(0.6666666666666666)
  b2 = T(0.8333333333333334)
  bbar1 = T(0.3333333333333333)
  bbar2 = T(0.4166666666666667)
  IRKN3ConstantCache(bconst1,bconst2,c1,a21,b1,b2,bbar1,bbar2)
end

function IRKN3ConstantCache(T::Type,T2::Type)
  bconst1 = T(3//2)
  bconst2 = T(-1//2)
  c1      = T2(1//2)
  a21     = T(1//8)
  b1      = T(2//3)
  b2      = T(5//6)
  bbar1   = T(1//3)
  bbar2   = T(5//12)
  IRKN3ConstantCache(bconst1,bconst2,c1,a21,b1,b2,bbar1,bbar2)
end

struct IRKN4ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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

Base.@pure function IRKN4ConstantCache{T<:CompiledFloats,T2<:CompiledFloats}(::Type{T},::Type{T2})
  bconst1 = T(1.5)
  bconst2 = T(-0.5)
  c1      = T2(0.25)
  c2      = T2(0.75)
  a21     = T(0.03125)
  # a31     = T(0)
  a32     = T(0.28125)
  b1      = T(1.0555555555555556)
  b2      = T(-0.16666666666666666)
  b3      = T(0.6111111111111112)
  bbar1   = T(-0.05555555555555555)
  bbar2   = T(0.2916666666666667)
  bbar3   = T(0.125)
  IRKN4ConstantCache(bconst1,bconst2,c1,c2,a21,a32,b1,b2,b3,bbar1,bbar2,bbar3)
end

function IRKN4ConstantCache(T::Type,T2::Type)
  bconst1 = T(3//2)
  bconst2 = T(-1//2)
  c1      = T2(1//4)
  c2      = T2(3//4)
  a21     = T(1//32)
  # a31     = T(0)
  a32     = T(9//32)
  b1      = T(19//18)
  b2      = T(-1//6)
  b3      = T(11//18)
  bbar1   = T(-1//18)
  bbar2   = T(7//24)
  bbar3   = T(1//8)
  IRKN4ConstantCache(bconst1,bconst2,c1,c2,a21,a32,b1,b2,b3,bbar1,bbar2,bbar3)
end

struct Nystrom5VelocityIndependentConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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

function Nystrom5VelocityIndependentConstantCache{T<:CompiledFloats,T2<:CompiledFloats}(::Type{T},::Type{T2})
  c1    = T2(0.2)
  c2    = T2(0.6666666666666666)
  # c3    = T2(1)
  a21   = T(0.02)
  a31   = T(-0.037037037037037035)
  a32   = T(0.25925925925925924)
  a41   = T(0.3)
  a42   = T(-0.05714285714285714)
  a43   = T(0.2571428571428571)
  bbar1 = T(0.041666666666666664)
  bbar2 = T(0.2976190476190476)
  bbar3 = T(0.16071428571428573)
  b1    = bbar1
  b2    = T(0.37202380952380953)
  b3    = T(0.48214285714285715)
  b4    = T(0.10416666666666667)
  Nystrom5VelocityIndependentConstantCache(c1, c2, a21, a31, a32, a41, a42, a43, bbar1, bbar2, bbar3, b1, b2, b3, b4)
end

function Nystrom5VelocityIndependentConstantCache(T::Type,T2::Type)
  c1    = T2(1//5)
  c2    = T2(2//3)
  # c3    = T2(1)
  a21   = T(1//50)
  a31   = T(-1//27)
  a32   = T(7//27)
  a41   = T(3//10)
  a42   = T(-2//35)
  a43   = T(9//35)
  bbar1 = T(14//336)
  bbar2 = T(100//336)
  bbar3 = T(54//336)
  b1    = bbar1
  b2    = T(125//336)
  b3    = T(162//336)
  b4    = T(35//336)
  Nystrom5VelocityIndependentConstantCache(c1, c2, a21, a31, a32, a41, a42, a43, bbar1, bbar2, bbar3, b1, b2, b3, b4)
end

struct DPRKN6ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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
end

Base.@pure function DPRKN6ConstantCache{T<:CompiledFloats,T2<:CompiledFloats}(::Type{T},::Type{T2})
  c1 = T2(0.12929590313670442)
  c2 = T2(0.25859180627340883)
  c3 = T2(0.67029708261548)
  c4 = T2(0.9)
  c5 = T2(1.0)
  a21 = T(0.008358715283968025)
  a31 = T(0.011144953711957367)
  a32 = T(0.022289907423914734)
  a41 = T(0.1454747428010918)
  a42 = T(-0.22986064052264749)
  a43 = T(0.3090349872029675)
  a51 = T(-0.20766826295078997)
  a52 = T(0.6863667842925143)
  a53 = T(-0.19954927787234925)
  a54 = T(0.12585075653062489)
  a61 = T(0.07811016144349478)
  a63 = T(0.2882917411897668)
  a64 = T(0.12242553717457041)
  a65 = T(0.011172560192168035)
  b1 = T(0.07811016144349478)
  b3 = T(0.2882917411897668)
  b4 = T(0.12242553717457041)
  b5 = T(0.011172560192168035)
  bp1 = T(0.07811016144349478)
  bp3 = T(0.3888434787059826)
  bp4 = T(0.3713207579288423)
  bp5 = T(0.11172560192168035)
  bp6 = T(0.05)
  btilde1 = T(-0.9807490989269235)
  btilde2 = T(2.406751371924452)
  btilde3 = T(-1.559600370364267)
  btilde4 = T(0.12242553717457041)
  btilde5 = T(0.011172560192168035)
  bptilde1 = T(0.023504273504273504)
  bptilde3 = T(-0.07242330719764424)
  bptilde4 = T(0.17543989844952962)
  bptilde5 = T(-0.2765208647561589)
  bptilde6 = T(0.15)
  DPRKN6ConstantCache(c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1, bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4, btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6)
end

function DPRKN6ConstantCache(T::Type,T2::Type)
  R        = sqrt(big(8581))
  c1       = T2((209-R)/900)
  c2       = T2((209-R)/450)
  c3       = T2((209+R)/450)
  c4       = T2(9//10)
  c5       = T2(1)
  a21      = T((26131-209R)/81_0000)
  a31      = T((26131-209R)/60_7500)
  a32      = T((26131-209R)/30_3750)
  a41      = T((980403512254+7781688431R)/116944_6992_1875)
  a42      = T(-(126288_4486208+153854_81287R)/116944_6992_1875)
  a43      = T((7166_233_891_441+786_945_632_99R)/46_777_879_687_500)
  a51      = T(-9(329260+3181R)/2704_0000)
  a52      = T(27(35129+3331R)/1352_0000)
  a53      = T(-27(554358343+31040327R)/46406048_0000)
  a54      = T(153(8555_257-67973R)/274592_0000)
  a61      = T(329//4212)
  # a62      = T(0)
  a63      = T((8411_9543+366_727R)/4096_22616)
  a64      = T((8411_9543-366_727R)/4096_22616)
  a65      = T(200//17901)
  b1       = T(329//4212)
  # b2       = T(0)
  b3       = a63
  b4       = a64
  b5       = T(200//17901)
  # b6       = T(0)
  bp1      = b1
  # bp2      = b2
  bp3      = T((389225579+96856R)/10_2405_6540)
  bp4      = T((389225579-96856R)/10_2405_6540)
  bp5      = T(2000//17901)
  bp6      = T(1//20)
  btilde1  = T(329//4212 - (2701+23R)/4563)
  btilde2  = T((9829+131R)/9126)
  btilde3  = T((8411_9543+366_727R)/4096_22616 - 5(1798+17R)/9126)
  btilde4  = b4
  btilde5  = b5
  # btilde6  = T(0)
  bptilde1 = T(329//4212 - 115//2106)
  # btildep2 = T(0)
  bptilde3 = T((389225579+96856R)/10_2405_6540 - (8411_9543+366_727R)/2560_14135)
  bptilde4 = T((389225579-96856R)/10_2405_6540 - (8411_9543-366_727R)/2560_14135)
  bptilde5 = T(2000//17901 - 6950//17901)
  bptilde6 = T(1//20 + 1//10)
  DPRKN6ConstantCache(c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1, bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4, btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6)
end

