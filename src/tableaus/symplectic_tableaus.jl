struct Symplectic2ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  a1::T
  a2::T
  b1::T
  b2::T
end

Base.@pure function PseudoVerletLeapfrogConstantCache(T,T2)
  a1 = T(1)
  a2 = T(0)
  b1 = T(1//2)
  b2 = T(1//2)
  Symplectic2ConstantCache{T,T2}(a1,a2,b1,b2)
end

Base.@pure function McAte2ConstantCache(T,T2)
  a2 = T(1 - (1/2)*sqrt(T(2)))
  a1 = T(1 - a2)
  b2 = T(1/(2*(1-a2)))
  b1 = T(1 - b2)
  Symplectic2ConstantCache{T,T2}(a1,a2,b1,b2)
end

Base.@pure function VerletLeapfrogConstantCache(T,T2)
  a1 = T(1//2)
  a2 = T(1//2)
  b1 = T(0)
  b2 = T(1)
  Symplectic2ConstantCache{T,T2}(a1,a2,b1,b2)
end

struct Symplectic3ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  a1::T
  a2::T
  a3::T
  b1::T
  b2::T
  b3::T
end

Base.@pure function Ruth3ConstantCache(T,T2)
  a1 = T(2//3)
  a2 = T(-2//3)
  a3 = T(1)
  b1 = T(7//24)
  b2 = T(3//4)
  b3 = T(-1//24)
  Symplectic3ConstantCache{T,T2}(a1,a2,a3,b1,b2,b3)
end

Base.@pure function McAte3ConstantCache(T,T2)
  a1 = T(0.9196615230173999)
  a2 = T(0.25/a1 - a1/2)
  a3 = T(1 - a1 - a2)
  b1 = T(a3)
  b2 = T(a2)
  b3 = T(a1)
  Symplectic3ConstantCache{T,T2}(a1,a2,a3,b1,b2,b3)
end

struct Symplectic4ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  a1::T
  a2::T
  a3::T
  a4::T
  b1::T
  b2::T
  b3::T
  b4::T
end

Base.@pure function CandyRoz4ConstantCache(T,T2)
  a1 = T((2 + T(2)^(1//3) + T(2)^(-1//3))/6)
  a2 = T((1 - T(2)^(1//3) - T(2)^(-1//3))/6)
  a3 = T(a2)
  a4 = T(a1)
  b1 = T(0)
  b2 = T((2 - T(2)^(1//3))^-1)
  b3 = T((1-T(2)^(2//3))^-1)
  b4 = T(b2)
  Symplectic4ConstantCache{T,T2}(a1,a2,a3,a4,b1,b2,b3,b4)
end

Base.@pure function McAte4ConstantCache(::Type{T},::Type{T2}) where {T<:CompiledFloats,T2<:CompiledFloats}
  a1 =T( 0.515352837431122936)
  a2 =T(-0.085782019412973646)
  a3 =T( 0.441583023616466524)
  a4 =T( 0.128846158365384185)
  b1 =T( 0.134496199277431089)
  b2 =T(-0.224819803079420806)
  b3 =T( 0.756320000515668291)
  b4 =T( 0.334003603286321425)
  Symplectic4ConstantCache{T,T2}(a1,a2,a3,a4,b1,b2,b3,b4)
end

Base.@pure function McAte4ConstantCache(T::Type,T2::Type)
  a1 =T(parse(BigFloat,"0.515352837431122936"))
  a2 =T(parse(BigFloat,"-0.085782019412973646"))
  a3 =T(parse(BigFloat," 0.441583023616466524"))
  a4 =T(parse(BigFloat," 0.128846158365384185"))
  b1 =T(parse(BigFloat," 0.134496199277431089"))
  b2 =T(parse(BigFloat,"-0.224819803079420806"))
  b3 =T(parse(BigFloat," 0.756320000515668291"))
  b4 =T(parse(BigFloat," 0.334003603286321425"))
  Symplectic4ConstantCache{T,T2}(a1,a2,a3,a4,b1,b2,b3,b4)
end

struct Symplectic45ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  a1::T
  a2::T
  a3::T
  a4::T
  a5::T
  b1::T
  b2::T
  b3::T
  b4::T
  b5::T
end

Base.@pure function CalvoSanz4ConstantCache(T,T2)
  a1= T(0.205177661542290)
  a2= T(0.403021281604210)
  a3=-T(0.12092087633891)
  a4= T(0.512721933192410)
  a5= T(0.0)
  b1= T(0.061758858135626)
  b2= T(0.33897802655364)
  b3= T(0.61479130717558)
  b4=-T(0.14054801465937)
  b5= T(0.12501982279453)
  Symplectic45ConstantCache{T,T2}(a1,a2,a3,a4,a5,b1,b2,b3,b4,b5)
end

# Broken
# http://epubs.siam.org/doi/pdf/10.1137/0916010
# On the numerical integration of ordinary differential equations by symmetric composition methods
Base.@pure function McAte42ConstantCache(::Type{T},::Type{T2}) where {T<:CompiledFloats,T2<:CompiledFloats}
  a1 =T(0.40518861839525227722)
  a2 =T(-0.28714404081652408900)
  a3 = 1 - 2a1 - 2a2
  a4 = a2
  a5 = a1
  b1 = T(-3//73)
  b2 = T(17//59)
  b3 = 1 - 2b1 - 2b2
  b4 = b2
  b5 = b1
  Symplectic45ConstantCache{T,T2}(a1,a2,a3,a4,a5,b1,b2,b3,b4,b5)
end

Base.@pure function McAte42ConstantCache(T::Type,T2::Type)
  a1 =T(parse(BigFloat,"0.40518861839525227722"))
  a2 =T(parse(BigFloat,"-0.28714404081652408900"))
  a3 = 1 - 2a1 - 2a2
  a4 = a2
  a5 = a1
  b1 = T(-3//73)
  b2 = T(17//59)
  b3 = 1 - 2b1 - 2b2
  b4 = b2
  b5 = b1
  Symplectic45ConstantCache{T,T2}(a1,a2,a3,a4,a5,b1,b2,b3,b4,b5)
end

struct Symplectic5ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  a1::T
  a2::T
  a3::T
  a4::T
  a5::T
  a6::T
  b1::T
  b2::T
  b3::T
  b4::T
  b5::T
  b6::T
end

Base.@pure function McAte5ConstantCache(::Type{T},::Type{T2}) where {T<:CompiledFloats,T2<:CompiledFloats}
  a1=T(0.339839625839110000)
  a2=T(-0.088601336903027329)
  a3=T(0.5858564768259621188)
  a4=T(-0.603039356536491888)
  a5=T(0.3235807965546976394)
  a6=T(0.4423637942197494587)
  b1=T(0.1193900292875672758)
  b2=T(0.6989273703824752308)
  b3=T(-0.1713123582716007754)
  b4=T(0.4012695022513534480)
  b5=T(0.0107050818482359840)
  b6=T(-0.0589796254980311632)
  Symplectic5ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6)
end

Base.@pure function McAte5ConstantCache(T::Type,T2::Type)
  a1=T(parse(BigFloat,"0.339839625839110000"))
  a2=T(parse(BigFloat,"-0.088601336903027329"))
  a3=T(parse(BigFloat,"0.5858564768259621188"))
  a4=T(parse(BigFloat,"-0.603039356536491888"))
  a5=T(parse(BigFloat,"0.3235807965546976394"))
  a6=T(parse(BigFloat,"0.4423637942197494587"))
  b1=T(parse(BigFloat,"0.1193900292875672758"))
  b2=T(parse(BigFloat,"0.6989273703824752308"))
  b3=T(parse(BigFloat,"-0.1713123582716007754"))
  b4=T(parse(BigFloat,"0.4012695022513534480"))
  b5=T(parse(BigFloat,"0.0107050818482359840"))
  b6=T(parse(BigFloat,"-0.0589796254980311632"))
  Symplectic5ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6)
end

struct Symplectic6ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  a1::T
  a2::T
  a3::T
  a4::T
  a5::T
  a6::T
  a7::T
  a8::T
  b1::T
  b2::T
  b3::T
  b4::T
  b5::T
  b6::T
  b7::T
  b8::T
end

Base.@pure function Yoshida6ConstantCache(T,T2)
  a1=T( 0.78451361047756)
  a2=T( 0.23557321335936)
  a3=T(-1.1776799841789)
  a4=T( 1.3151863206839)
  a5=T( a3)
  a6=T( a2)
  a7=T( a1)
  a8=T( 0.0)
  b1= a1/2
  b2=T((a1+a2)/2)
  b3=T((a2+a3)/2)
  b4=T((a3+a4)/2)
  b5=T( b4)
  b6=T( b3)
  b7=T( b2)
  b8=T( b1)
  Symplectic6ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8)
end

struct Symplectic62ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  a1::T
  a2::T
  a3::T
  a4::T
  a5::T
  a6::T
  a7::T
  a8::T
  a9::T
  a10::T
  b1::T
  b2::T
  b3::T
  b4::T
  b5::T
  b6::T
  b7::T
  b8::T
  b9::T
  b10::T
end

Base.@pure function KahanLi6ConstantCache(::Type{T},::Type{T2}) where {T<:CompiledFloats,T2<:CompiledFloats}
  a1=T(0.39216144400731413927925056)
  a2=T(0.33259913678935943859974864)
  a3=T(-0.70624617255763935980996482)
  a4=T(0.08221359629355080023149045)
  a5=T(0.79854399093482996339895035)
  a6= a4
  a7= a3
  a8= a2
  a9= a1
  a10 = T(0)
  b1= T(a1/2)
  b2=T((a1+a2)/2)
  b3=T((a2+a3)/2)
  b4=T((a3+a4)/2)
  b5=T((a4+a5)/2)
  b6=b5
  b7=b4
  b8=b3
  b9=b2
  b10=b1
  Symplectic62ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
end

Base.@pure function KahanLi6ConstantCache(T::Type,T2::Type)
  a1=T(parse(BigFloat,"0.39216144400731413927925056"))
  a2=T(parse(BigFloat,"0.33259913678935943859974864"))
  a3=T(parse(BigFloat,"-0.70624617255763935980996482"))
  a4=T(parse(BigFloat,"0.08221359629355080023149045"))
  a5=T(parse(BigFloat,"0.79854399093482996339895035"))
  a6= a4
  a7= a3
  a8= a2
  a9= a1
  a10 = T(0)
  b1= T(a1/2)
  b2=T((a1+a2)/2)
  b3=T((a2+a3)/2)
  b4=T((a3+a4)/2)
  b5=T((a4+a5)/2)
  b6=b5
  b7=b4
  b8=b3
  b9=b2
  b10=b1
  Symplectic62ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
end

struct McAte8ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  a1::T
  a2::T
  a3::T
  a4::T
  a5::T
  a6::T
  a7::T
  a8::T
  a9::T
  a10::T
  a11::T
  a12::T
  a13::T
  a14::T
  a15::T
  a16::T
  b1::T
  b2::T
  b3::T
  b4::T
  b5::T
  b6::T
  b7::T
  b8::T
  b9::T
  b10::T
  b11::T
  b12::T
  b13::T
  b14::T
  b15::T
  b16::T
end

Base.@pure function McAte8ConstantCache(::Type{T},::Type{T2}) where {T<:CompiledFloats,T2<:CompiledFloats}
  a1=T(0.74167036435061295344822780)
  a2=T(-0.40910082580003159399730010)
  a3=T(0.19075471029623837995387626)
  a4=T(-0.57386247111608226665638773)
  a5=T(0.29906418130365592384446354)
  a6=T(0.33462491824529818378495798)
  a7=T(0.31529309239676659663205666)
  a8=T(-0.79688793935291635401978884)
  a9=a7
  a10=a6
  a11=a5
  a12=a4
  a13=a3
  a14=a2
  a15=a1
  a16=T(0)
  b1= T(a1/2)
  b2=T((a1+a2)/2)
  b3=T((a2+a3)/2)
  b4=T((a3+a4)/2)
  b5=T((a4+a5)/2)
  b6=T((a5+a6)/2)
  b7=T((a6+a7)/2)
  b8=T((a7+a8)/2)
  b9=b8
  b10=b7
  b11=b6
  b12=b5
  b13=b4
  b14=b3
  b15=b2
  b16=b1
  McAte8ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,
                            b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16)
end

Base.@pure function McAte8ConstantCache(T::Type,T2::Type)
  a1=T(parse(BigFloat,"0.74167036435061295344822780"))
  a2=T(parse(BigFloat,"-0.40910082580003159399730010"))
  a3=T(parse(BigFloat,"0.19075471029623837995387626"))
  a4=T(parse(BigFloat,"-0.57386247111608226665638773"))
  a5=T(parse(BigFloat,"0.29906418130365592384446354"))
  a6=T(parse(BigFloat,"0.33462491824529818378495798"))
  a7=T(parse(BigFloat,"0.31529309239676659663205666"))
  a8=T(parse(BigFloat,"-0.79688793935291635401978884"))
  a9=a7
  a10=a6
  a11=a5
  a12=a4
  a13=a3
  a14=a2
  a15=a1
  a16=T(0)
  b1= T(a1/2)
  b2=T((a1+a2)/2)
  b3=T((a2+a3)/2)
  b4=T((a3+a4)/2)
  b5=T((a4+a5)/2)
  b6=T((a5+a6)/2)
  b7=T((a6+a7)/2)
  b8=T((a7+a8)/2)
  b9=b8
  b10=b7
  b11=b6
  b12=b5
  b13=b4
  b14=b3
  b15=b2
  b16=b1
  McAte8ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,
                            b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16)
end

struct KahanLi8ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  a1::T
  a2::T
  a3::T
  a4::T
  a5::T
  a6::T
  a7::T
  a8::T
  a9::T
  a10::T
  a11::T
  a12::T
  a13::T
  a14::T
  a15::T
  a16::T
  a17::T
  a18::T
  b1::T
  b2::T
  b3::T
  b4::T
  b5::T
  b6::T
  b7::T
  b8::T
  b9::T
  b10::T
  b11::T
  b12::T
  b13::T
  b14::T
  b15::T
  b16::T
  b17::T
  b18::T
end

Base.@pure function KahanLi8ConstantCache(::Type{T},::Type{T2}) where {T<:CompiledFloats,T2<:CompiledFloats}
  a1=T(0.13020248308889008087881763)
  a2=T(0.56116298177510838456196441)
  a3=T(-0.38947496264484728640807860)
  a4=T(0.15884190655515560089621075)
  a5=T(-0.39590389413323757733623154)
  a6=T(0.18453964097831570709183254)
  a7=T(0.25837438768632204729397911)
  a8=T(0.29501172360931029887096624)
  a9=T(-0.60550853383003451169892108)
  a10=a8
  a11=a7
  a12=a6
  a13=a5
  a14=a4
  a15=a3
  a16=a2
  a17=a1
  a18=T(0)
  b1= T(a1/2)
  b2=T((a1+a2)/2)
  b3=T((a2+a3)/2)
  b4=T((a3+a4)/2)
  b5=T((a4+a5)/2)
  b6=T((a5+a6)/2)
  b7=T((a6+a7)/2)
  b8=T((a7+a8)/2)
  b9=T((a8+a9)/2)
  b10=b9
  b11=b8
  b12=b7
  b13=b6
  b14=b5
  b15=b4
  b16=b3
  b17=b2
  b18=b1
  KahanLi8ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
                              b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18)
end

Base.@pure function KahanLi8ConstantCache(T::Type,T2::Type)
  a1=T(parse(BigFlaot,"0.13020248308889008087881763"))
  a2=T(parse(BigFlaot,"0.56116298177510838456196441"))
  a3=T(parse(BigFlaot,"-0.38947496264484728640807860"))
  a4=T(parse(BigFlaot,"0.15884190655515560089621075"))
  a5=T(parse(BigFlaot,"-0.39590389413323757733623154"))
  a6=T(parse(BigFlaot,"0.18453964097831570709183254"))
  a7=T(parse(BigFlaot,"0.25837438768632204729397911"))
  a8=T(parse(BigFlaot,"0.29501172360931029887096624"))
  a9=T(parse(BigFlaot,"-0.60550853383003451169892108"))
  a10=a8
  a11=a7
  a12=a6
  a13=a5
  a14=a4
  a15=a3
  a16=a2
  a17=a1
  a18=T(0)
  b1= T(a1/2)
  b2=T((a1+a2)/2)
  b3=T((a2+a3)/2)
  b4=T((a3+a4)/2)
  b5=T((a4+a5)/2)
  b6=T((a5+a6)/2)
  b7=T((a6+a7)/2)
  b8=T((a7+a8)/2)
  b9=T((a8+a9)/2)
  b10=b9
  b11=b8
  b12=b7
  b13=b6
  b14=b5
  b15=b4
  b16=b3
  b17=b2
  b18=b1
  KahanLi8ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
                              b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18)
end

struct SofSpa10ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  a1::T
  a2::T
  a3::T
  a4::T
  a5::T
  a6::T
  a7::T
  a8::T
  a9::T
  a10::T
  a11::T
  a12::T
  a13::T
  a14::T
  a15::T
  a16::T
  a17::T
  a18::T
  a19::T
  a20::T
  a21::T
  a22::T
  a23::T
  a24::T
  a25::T
  a26::T
  a27::T
  a28::T
  a29::T
  a30::T
  a31::T
  a32::T
  a33::T
  a34::T
  a35::T
  a36::T
  b1::T
  b2::T
  b3::T
  b4::T
  b5::T
  b6::T
  b7::T
  b8::T
  b9::T
  b10::T
  b11::T
  b12::T
  b13::T
  b14::T
  b15::T
  b16::T
  b17::T
  b18::T
  b19::T
  b20::T
  b21::T
  b22::T
  b23::T
  b24::T
  b25::T
  b26::T
  b27::T
  b28::T
  b29::T
  b30::T
  b31::T
  b32::T
  b33::T
  b34::T
  b35::T
  b36::T
end

Base.@pure function SofSpa10ConstantCache(::Type{T},::Type{T2}) where {T<:CompiledFloats,T2<:CompiledFloats}
  a1=T(0.07879572252168641926390768)
  a2=T(0.31309610341510852776481247)
  a3=T(0.02791838323507806610952027)
  a4=T(-0.22959284159390709415121340)
  a5=T(0.13096206107716486317465686)
  a6=T(-0.26973340565451071434460973)
  a7=T(0.07497334315589143566613711)
  a8=T(0.11199342399981020488957508)
  a9=T(0.36613344954622675119314812)
  a10=T(-0.39910563013603589787862981)
  a11=T(0.10308739852747107731580277)
  a12=T(0.41143087395589023782070412)
  a13=T(-0.00486636058313526176219566)
  a14=T(-0.39203335370863990644808194)
  a15=T(0.05194250296244964703718290)
  a16=T(0.05066509075992449633587434)
  a17=T(0.04967437063972987905456880)
  a18=T(0.04931773575959453791768001)
  a19=a17
  a20=a16
  a21=a15
  a22=a14
  a23=a13
  a24=a12
  a25=a11
  a26=a10
  a27=a9
  a28=a8
  a29=a7
  a30=a6
  a31=a5
  a32=a4
  a33=a3
  a34=a2
  a35=a1
  a36=T(0)
  b1= T(a1/2)
  b2=T((a1+a2)/2)
  b3=T((a2+a3)/2)
  b4=T((a3+a4)/2)
  b5=T((a4+a5)/2)
  b6=T((a5+a6)/2)
  b7=T((a6+a7)/2)
  b8=T((a7+a8)/2)
  b9=T((a8+a9)/2)
  b10=T((a9+a10)/2)
  b11=T((a10+a11)/2)
  b12=T((a11+a12)/2)
  b13=T((a12+a13)/2)
  b14=T((a13+a14)/2)
  b15=T((a14+a15)/2)
  b16=T((a15+a16)/2)
  b17=T((a16+a17)/2)
  b18=T((a17+a18)/2)
  b19=b18
  b20=b17
  b21=b16
  b22=b15
  b23=b14
  b24=b13
  b25=b12
  b26=b11
  b27=b10
  b28=b9
  b29=b8
  b30=b7
  b31=b6
  b32=b5
  b33=b4
  b34=b3
  b35=b2
  b36=b1
  SofSpa10ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
                              a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,
                              a35,a36,
                              b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,
                              b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,
                              b35,b36)
end

Base.@pure function SofSpa10ConstantCache(T::Type,T2::Type)
  a1=T(parse(BigFloat,"0.07879572252168641926390768"))
  a2=T(parse(BigFloat,"0.31309610341510852776481247"))
  a3=T(parse(BigFloat,"0.02791838323507806610952027"))
  a4=T(parse(BigFloat,"-0.22959284159390709415121340"))
  a5=T(parse(BigFloat,"0.13096206107716486317465686"))
  a6=T(parse(BigFloat,"-0.26973340565451071434460973"))
  a7=T(parse(BigFloat,"0.07497334315589143566613711"))
  a8=T(parse(BigFloat,"0.11199342399981020488957508"))
  a9=T(parse(BigFloat,"0.36613344954622675119314812"))
  a10=T(parse(BigFloat,"-0.39910563013603589787862981"))
  a11=T(parse(BigFloat,"0.10308739852747107731580277"))
  a12=T(parse(BigFloat,"0.41143087395589023782070412"))
  a13=T(parse(BigFloat,"-0.00486636058313526176219566"))
  a14=T(parse(BigFloat,"-0.39203335370863990644808194"))
  a15=T(parse(BigFloat,"0.05194250296244964703718290"))
  a16=T(parse(BigFloat,"0.05066509075992449633587434"))
  a17=T(parse(BigFloat,"0.04967437063972987905456880"))
  a18=T(parse(BigFloat,"0.04931773575959453791768001"))
  a19=a17
  a20=a16
  a21=a15
  a22=a14
  a23=a13
  a24=a12
  a25=a11
  a26=a10
  a27=a9
  a28=a8
  a29=a7
  a30=a6
  a31=a5
  a32=a4
  a33=a3
  a34=a2
  a35=a1
  a36=T(0)
  b1= T(a1/2)
  b2=T((a1+a2)/2)
  b3=T((a2+a3)/2)
  b4=T((a3+a4)/2)
  b5=T((a4+a5)/2)
  b6=T((a5+a6)/2)
  b7=T((a6+a7)/2)
  b8=T((a7+a8)/2)
  b9=T((a8+a9)/2)
  b10=T((a9+a10)/2)
  b11=T((a10+a11)/2)
  b12=T((a11+a12)/2)
  b13=T((a12+a13)/2)
  b14=T((a13+a14)/2)
  b15=T((a14+a15)/2)
  b16=T((a15+a16)/2)
  b17=T((a16+a17)/2)
  b18=T((a17+a18)/2)
  b19=b18
  b20=b17
  b21=b16
  b22=b15
  b23=b14
  b24=b13
  b25=b12
  b26=b11
  b27=b10
  b28=b9
  b29=b8
  b30=b7
  b31=b6
  b32=b5
  b33=b4
  b34=b3
  b35=b2
  b36=b1
  SofSpa10ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
                              a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,
                              a35,a36,
                              b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,
                              b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,
                              b35,b36)
end
