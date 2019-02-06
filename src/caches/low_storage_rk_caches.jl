
@cache struct CarpenterKennedy2N54Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct CarpenterKennedy2N54ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  A2end::SVector{4,T} # A1 is always zero
  B1::T
  B2end::SVector{4,T}
  c2end::SVector{4,T2} # c1 is always zero

  function CarpenterKennedy2N54ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    A2 = T(-567301805773/1357537059087)
    A3 = T(-2404267990393/2016746695238)
    A4 = T(-3550918686646/2091501179385)
    A5 = T(-1275806237668/842570457699)
    A2end = SVector{4,T}(A2, A3, A4, A5)

    B1 = T(1432997174477/9575080441755)
    B2 = T(5161836677717/13612068292357)
    B3 = T(1720146321549/2090206949498)
    B4 = T(3134564353537/4481467310338)
    B5 = T(2277821191437/14882151754819)
    B2end = SVector{4,T}(B2, B3, B4, B5)

    c2 = T2(1432997174477/9575080441755)
    c3 = T2(2526269341429/6820363962896)
    c4 = T2(2006345519317/3224310063776)
    c5 = T2(2802321613138/2924317926251)
    c2end = SVector{4,T2}(c2, c3, c4, c5)

    new{T,T2}(A2end, B1, B2end, c2end)
  end
end

function alg_cache(alg::CarpenterKennedy2N54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = CarpenterKennedy2N54ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  CarpenterKennedy2N54Cache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::CarpenterKennedy2N54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  CarpenterKennedy2N54ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end
