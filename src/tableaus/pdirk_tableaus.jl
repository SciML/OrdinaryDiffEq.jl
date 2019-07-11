struct PDIRK44ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  γs::SVector{T2}
  cs::SVector{T2}
  α1::SVector{T}
  α2::SVector{T}
end

