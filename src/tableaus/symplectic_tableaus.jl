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
