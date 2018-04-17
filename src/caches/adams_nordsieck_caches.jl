mutable struct AN5ConstantCache{histType,lType,dtType} <: OrdinaryDiffEqConstantCache
  hist1::histType
  hist2::histType
  l::MVector{6,lType}
  m::MVector{6,lType}
  tq::lType
  tau::MVector{6, dtType}
end

function AN5ConstantCache(rate_prototype, tTypeNoUnits, dt)
  siz = size(rate_prototype)
  typ = Base.promote_op(*, eltype(rate_prototype), tTypeNoUnits)
  hist1 = zeros(typ, siz); hist2 = zeros(hist1)
  l = zeros(tTypeNoUnits, MVector{6}); m = zeros(l)
  tq = zero(tTypeNoUnits)
  tau = zeros(dt, MVector{6})
  AN5ConstantCache(hist1,hist2,l,m,tq, tau)
end

function alg_cache(alg::AN5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  AN5ConstantCache(rate_prototype, tTypeNoUnits, dt)
end

mutable struct AN5Cache{uType,rateType,histType,lType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  tmp::uType
  tab::AN5ConstantCache{histType,lType}
end

u_cache(c::AN5Cache) = ()
du_cache(c::AN5Cache) = (c.fsalfirst,c.hist1,c.hist2)

function alg_cache(alg::AN5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  fsalfirst = zeros(rate_prototype)
  tmp = similar(u)
  tab = AN5ConstantCache(rate_prototype, tTypeNoUnits, dt)
  AN5Cache(u,uprev,fsalfirst,tmp,tab)
end
