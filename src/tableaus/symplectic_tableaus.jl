immutable Symplectic2ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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

immutable Symplectic3ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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

immutable Symplectic4ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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
  a1 = T((2 + 2^(1//3) + 2^(-1//3))/6)
  a2 = T((1 - 2^(1//3) - 2^(-1//3))/6)
  a3 = T(a2)
  a4 = T(a1)
  b1 = T(0)
  b2 = T((2 - 2^(1//3))^-1)
  b3 = T((1-2^(2//3))^-1)
  b4 = T(b2)
  Symplectic4ConstantCache{T,T2}(a1,a2,a3,a4,b1,b2,b3,b4)
end

Base.@pure function McAte4ConstantCache(T,T2)
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

immutable Symplectic45ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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

immutable Symplectic5ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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

Base.@pure function McAte5ConstantCache(T,T2)
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

immutable Symplectic6ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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
  b1=T( 0.39225680523878)
  b2=T( 0.51004341191846)
  b3=T(-0.47105338540976)
  b4=T( 0.068753168252520)
  b5=T( b4)
  b6=T( b3)
  b7=T( b2)
  b8=T( b1)
  Symplectic6ConstantCache{T,T2}(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8)
end
