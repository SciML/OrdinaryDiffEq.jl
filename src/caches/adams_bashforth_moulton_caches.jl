@cache mutable struct AB3Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  ralk2::rateType
  k::rateType
  tmp::uType
  step::Int
end

@cache mutable struct AB3ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  step::Int
end

function alg_cache(alg::AB3,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  ralk2 = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = zero(u)
  AB3Cache(u,uprev,fsalfirst,k2,k3,ralk2,k,tmp,1)
end

function alg_cache(alg::AB3,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  k2 = rate_prototype
  k3 = rate_prototype
  AB3ConstantCache(k2,k3,1)
end

@cache mutable struct ABM32Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  ralk2::rateType
  k::rateType
  tmp::uType
  step::Int
end

@cache mutable struct ABM32ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  step::Int
end

function alg_cache(alg::ABM32,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  ralk2 = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = zero(u)
  ABM32Cache(u,uprev,fsalfirst,k2,k3,ralk2,k,tmp,1)
end

function alg_cache(alg::ABM32,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  k2 = rate_prototype
  k3 = rate_prototype
  ABM32ConstantCache(k2,k3,1)
end

@cache mutable struct AB4Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  ralk2::rateType
  k::rateType
  tmp::uType
  t2::rateType
  t3::rateType
  t4::rateType
  step::Int
end

@cache mutable struct AB4ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  k4::rateType
  step::Int
end

function alg_cache(alg::AB4,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  ralk2 = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = zero(u)
  t2 = zero(rate_prototype)
  t3 = zero(rate_prototype)
  t4 = zero(rate_prototype)
  AB4Cache(u,uprev,fsalfirst,k2,k3,k4,ralk2,k,tmp,t2,t3,t4,1)
end

function alg_cache(alg::AB4,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  k2 = rate_prototype
  k3 = rate_prototype
  k4 = rate_prototype
  AB4ConstantCache(k2,k3,k4,1)
end

@cache mutable struct ABM43Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  ralk2::rateType
  k::rateType
  tmp::uType
  t2::rateType
  t3::rateType
  t4::rateType
  t5::rateType
  t6::rateType
  t7::rateType
  step::Int
end

@cache mutable struct ABM43ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  k4::rateType
  step::Int
end

function alg_cache(alg::ABM43,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  ralk2 = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = zero(u)
  t2 = zero(rate_prototype)
  t3 = zero(rate_prototype)
  t4 = zero(rate_prototype)
  t5 = zero(rate_prototype)
  t6 = zero(rate_prototype)
  t7 = zero(rate_prototype)
  ABM43Cache(u,uprev,fsalfirst,k2,k3,k4,ralk2,k,tmp,t2,t3,t4,t5,t6,t7,1)
end

function alg_cache(alg::ABM43,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  k2 = rate_prototype
  k3 = rate_prototype
  k4 = rate_prototype
  ABM43ConstantCache(k2,k3,k4,1)
end

@cache mutable struct AB5Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k::rateType
  tmp::uType
  t2::rateType
  t3::rateType
  t4::rateType
  step::Int
end

@cache mutable struct AB5ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  step::Int
end

function alg_cache(alg::AB5,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = zero(u)
  t2 = zero(rate_prototype)
  t3 = zero(rate_prototype)
  t4 = zero(rate_prototype)
  AB5Cache(u,uprev,fsalfirst,k2,k3,k4,k5,k,tmp,t2,t3,t4,1)
end

function alg_cache(alg::AB5,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  k2 = rate_prototype
  k3 = rate_prototype
  k4 = rate_prototype
  k5 = rate_prototype
  AB5ConstantCache(k2,k3,k4,k5,1)
end

@cache mutable struct ABM54Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k::rateType
  tmp::uType
  t2::rateType
  t3::rateType
  t4::rateType
  t5::rateType
  t6::rateType
  t7::rateType
  t8::rateType
  step::Int
end

@cache mutable struct ABM54ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  step::Int
end

function alg_cache(alg::ABM54,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = zero(u)
  t2 = zero(rate_prototype)
  t3 = zero(rate_prototype)
  t4 = zero(rate_prototype)
  t5 = zero(rate_prototype)
  t6 = zero(rate_prototype)
  t7 = zero(rate_prototype)
  t8 = zero(rate_prototype)
  ABM54Cache(u,uprev,fsalfirst,k2,k3,k4,k5,k,tmp,t2,t3,t4,t5,t6,t7,t8,1)
end

function alg_cache(alg::ABM54,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  k2 = rate_prototype
  k3 = rate_prototype
  k4 = rate_prototype
  k5 = rate_prototype
  ABM54ConstantCache(k2,k3,k4,k5,1)
end

@cache mutable struct VCAB3ConstantCache{TabType,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::rArrayType
  ϕstar_nm1::rArrayType
  ϕstar_n::rArrayType
  β::tArrayType
  order::Int
  tab::TabType
  step::Int
end

@cache mutable struct VCAB3Cache{uType,rateType,TabType,bs3Type,tArrayType,cArrayType,uNoUnitsType,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  bs3cache::bs3Type
  k4::rateType
  ϕstar_nm1::coefType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::coefType
  ϕstar_n::coefType
  β::tArrayType
  order::Int
  atmp::uNoUnitsType
  tmp::uType
  utilde::uType
  tab::TabType
  step::Int
end

function alg_cache(alg::VCAB3,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  dts = fill(zero(dt),3)
  c = fill(zero(t), 3, 3)
  g = fill(zero(t), 3)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 3)
  for i in 1:3
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(t),3)
  order = 3
  tab = BS3ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  VCAB3ConstantCache(dts,c,g,ϕ_n,ϕstar_nm1,ϕstar_n,β,order,tab,1)
end

function alg_cache(alg::VCAB3,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = BS3ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  bk1 = zero(rate_prototype)
  bk2 = zero(rate_prototype)
  bk3 = zero(rate_prototype)
  bk4 = zero(rate_prototype)
  butilde = zero(u)
  batmp = similar(u,uEltypeNoUnits)
  btmp = zero(u)
  bs3cache = BS3Cache(u,uprev,bk1,bk2,bk3,bk4,butilde,btmp,batmp,tab)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(dt),3)
  c = fill(zero(t), 3, 3)
  g = fill(zero(t), 3)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 3)
  for i in 1:3
    ϕ_n[i] = zero(rate_prototype)
    ϕstar_nm1[i] = zero(rate_prototype)
    ϕstar_n[i] = zero(rate_prototype)
  end
  β = fill(zero(t),3)
  order = 3
  atmp = similar(u,uEltypeNoUnits)
  tmp = zero(u)
  utilde = zero(u)
  VCAB3Cache(u,uprev,fsalfirst,bs3cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,atmp,tmp,utilde,tab,1)
end

@cache mutable struct VCAB4ConstantCache{rk4constcache,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
  ϕstar_nm1::rArrayType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::rArrayType
  ϕstar_n::rArrayType
  β::tArrayType
  order::Int
  rk4constcache::rk4constcache
  step::Int
end

@cache mutable struct VCAB4Cache{uType,rateType,rk4cacheType,tArrayType,cArrayType,uNoUnitsType,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  rk4cache::rk4cacheType
  k4::rateType
  ϕstar_nm1::coefType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::coefType
  ϕstar_n::coefType
  β::tArrayType
  order::Int
  atmp::uNoUnitsType
  tmp::uType
  utilde::uType
  step::Int
end

function alg_cache(alg::VCAB4,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  dts = fill(zero(dt),4)
  c = fill(zero(t), 4, 4)
  g = fill(zero(t), 4)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 4)
  for i in 1:4
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(t),4)
  order = 4
  rk4constcache = RK4ConstantCache()
  VCAB4ConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,rk4constcache,1)
end

function alg_cache(alg::VCAB4,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  rk1 = zero(rate_prototype)
  rk2 = zero(rate_prototype)
  rk3 = zero(rate_prototype)
  rk4 = zero(rate_prototype)
  rk  = zero(rate_prototype)
  rtmp = zero(u); ratmp = similar(u, uEltypeNoUnits)
  rk4cache = RK4Cache(u,uprev,rk1,rk2,rk3,rk4,rk,rtmp,ratmp)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(dt),4)
  c = fill(zero(t), 4, 4)
  g = fill(zero(t), 4)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 4)
  for i in 1:4
    ϕ_n[i] = zero(rate_prototype)
    ϕstar_nm1[i] = zero(rate_prototype)
    ϕstar_n[i] = zero(rate_prototype)
  end
  β = fill(zero(t),4)
  order = 4
  atmp = similar(u,uEltypeNoUnits)
  tmp = zero(u)
  utilde = zero(u)
  VCAB4Cache(u,uprev,fsalfirst,rk4cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,atmp,tmp,utilde,1)
end

# VCAB5

@cache mutable struct VCAB5ConstantCache{rk4constcache,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
  ϕstar_nm1::rArrayType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::rArrayType
  ϕstar_n::rArrayType
  β::tArrayType
  order::Int
  rk4constcache::rk4constcache
  step::Int
end

@cache mutable struct VCAB5Cache{uType,rateType,rk4cacheType,tArrayType,cArrayType,uNoUnitsType,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  rk4cache::rk4cacheType
  k4::rateType
  ϕstar_nm1::coefType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::coefType
  ϕstar_n::coefType
  β::tArrayType
  order::Int
  atmp::uNoUnitsType
  tmp::uType
  utilde::uType
  step::Int
end

function alg_cache(alg::VCAB5,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  dts = fill(zero(dt),5)
  c = fill(zero(t), 5, 5)
  g = fill(zero(t), 5)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 5)
  for i in 1:5
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(t),5)
  order = 5
  rk4constcache = RK4ConstantCache()
  VCAB5ConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,rk4constcache,1)
end

function alg_cache(alg::VCAB5,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  rk1 = zero(rate_prototype)
  rk2 = zero(rate_prototype)
  rk3 = zero(rate_prototype)
  rk4 = zero(rate_prototype)
  rk  = zero(rate_prototype)
  rtmp = zero(u); ratmp = similar(u, uEltypeNoUnits)
  rk4cache = RK4Cache(u,uprev,rk1,rk2,rk3,rk4,rk,rtmp,ratmp)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(dt),5)
  c = fill(zero(t), 5, 5)
  g = fill(zero(t), 5)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 5)
  for i in 1:5
    ϕ_n[i] = zero(rate_prototype)
    ϕstar_nm1[i] = zero(rate_prototype)
    ϕstar_n[i] = zero(rate_prototype)
  end
  β = fill(zero(t),5)
  order = 5
  atmp = similar(u,uEltypeNoUnits)
  tmp = zero(u)
  utilde = zero(u)
  VCAB5Cache(u,uprev,fsalfirst,rk4cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,atmp,tmp,utilde,1)
end

# VCABM3

@cache mutable struct VCABM3ConstantCache{TabType,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::rArrayType
  ϕ_np1::rArrayType
  ϕstar_nm1::rArrayType
  ϕstar_n::rArrayType
  β::tArrayType
  order::Int
  tab::TabType
  step::Int
end

@cache mutable struct VCABM3Cache{uType,rateType,TabType,bs3Type,tArrayType,cArrayType,uNoUnitsType,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  bs3cache::bs3Type
  k4::rateType
  ϕstar_nm1::coefType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::coefType
  ϕ_np1::coefType
  ϕstar_n::coefType
  β::tArrayType
  order::Int
  atmp::uNoUnitsType
  tmp::uType
  utilde::uType
  tab::TabType
  step::Int
end

function alg_cache(alg::VCABM3,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  dts = fill(zero(dt),3)
  c = fill(zero(t), 4, 4)
  g = fill(zero(t), 4)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 3)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 4)
  for i in 1:3
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(t),3)
  order = 3
  tab = BS3ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  VCABM3ConstantCache(dts,c,g,ϕ_n,ϕ_np1,ϕstar_nm1,ϕstar_n,β,order,tab,1)
end

function alg_cache(alg::VCABM3,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = BS3ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  bk1 = zero(rate_prototype)
  bk2 = zero(rate_prototype)
  bk3 = zero(rate_prototype)
  bk4 = zero(rate_prototype)
  butilde = zero(u)
  batmp = similar(u,uEltypeNoUnits)
  btmp = zero(u)
  bs3cache = BS3Cache(u,uprev,bk1,bk2,bk3,bk4,butilde,btmp,batmp,tab)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(dt),3)
  c = fill(zero(t), 4, 4)
  g = fill(zero(t), 4)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 3)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 4)
  for i in 1:3
    ϕ_n[i] = zero(rate_prototype)
    ϕstar_nm1[i] = zero(rate_prototype)
    ϕstar_n[i] = zero(rate_prototype)
  end
  for i in 1:4
    ϕ_np1[i] = zero(rate_prototype)
  end
  β = fill(zero(t),3)
  order = 3
  atmp = similar(u,uEltypeNoUnits)
  tmp = zero(u)
  utilde = zero(u)
  VCABM3Cache(u,uprev,fsalfirst,bs3cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,atmp,tmp,utilde,tab,1)
end

# VCABM4

@cache mutable struct VCABM4ConstantCache{rk4constcache,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
  ϕstar_nm1::rArrayType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::rArrayType
  ϕ_np1::rArrayType
  ϕstar_n::rArrayType
  β::tArrayType
  order::Int
  rk4constcache::rk4constcache
  step::Int
end

@cache mutable struct VCABM4Cache{uType,rateType,rk4cacheType,tArrayType,cArrayType,uNoUnitsType,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  rk4cache::rk4cacheType
  k4::rateType
  ϕstar_nm1::coefType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::coefType
  ϕ_np1::coefType
  ϕstar_n::coefType
  β::tArrayType
  order::Int
  atmp::uNoUnitsType
  tmp::uType
  utilde::uType
  step::Int
end

function alg_cache(alg::VCABM4,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  dts = fill(zero(dt),4)
  c = fill(zero(t), 5, 5)
  g = fill(zero(t), 5)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 4)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 5)
  for i in 1:4
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(t),4)
  order = 4
  rk4constcache = RK4ConstantCache()
  VCABM4ConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,rk4constcache,1)
end

function alg_cache(alg::VCABM4,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  rk1 = zero(rate_prototype)
  rk2 = zero(rate_prototype)
  rk3 = zero(rate_prototype)
  rk4 = zero(rate_prototype)
  rk  = zero(rate_prototype)
  rtmp = zero(u); ratmp = similar(u, uEltypeNoUnits)
  rk4cache = RK4Cache(u,uprev,rk1,rk2,rk3,rk4,rk,rtmp,ratmp)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(dt),4)
  c = fill(zero(t), 5, 5)
  g = fill(zero(t), 5)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 4)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 5)
  for i in 1:4
    ϕ_n[i] = zero(rate_prototype)
    ϕstar_nm1[i] = zero(rate_prototype)
    ϕstar_n[i] = zero(rate_prototype)
  end
  for i in 1:5
    ϕ_np1[i] = zero(rate_prototype)
  end
  β = fill(zero(t),4)
  order = 4
  atmp = similar(u,uEltypeNoUnits)
  tmp = zero(u)
  utilde = zero(u)
  VCABM4Cache(u,uprev,fsalfirst,rk4cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,atmp,tmp,utilde,1)
end

# VCABM5

@cache mutable struct VCABM5ConstantCache{rk4constcache,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
  ϕstar_nm1::rArrayType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::rArrayType
  ϕ_np1::rArrayType
  ϕstar_n::rArrayType
  β::tArrayType
  order::Int
  rk4constcache::rk4constcache
  step::Int
end

@cache mutable struct VCABM5Cache{uType,rateType,rk4cacheType,tArrayType,cArrayType,uNoUnitsType,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  rk4cache::rk4cacheType
  k4::rateType
  ϕstar_nm1::coefType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::coefType
  ϕ_np1::coefType
  ϕstar_n::coefType
  β::tArrayType
  order::Int
  atmp::uNoUnitsType
  tmp::uType
  utilde::uType
  step::Int
end

function alg_cache(alg::VCABM5,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  dts = fill(zero(t),5)
  c = fill(zero(t), 6, 6)
  g = fill(zero(t), 6)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 5)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 6)
  for i in 1:5
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(t),5)
  order = 5
  rk4constcache = RK4ConstantCache()
  VCABM5ConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,rk4constcache,1)
end

function alg_cache(alg::VCABM5,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  rk1 = zero(rate_prototype)
  rk2 = zero(rate_prototype)
  rk3 = zero(rate_prototype)
  rk4 = zero(rate_prototype)
  rk  = zero(rate_prototype)
  rtmp = zero(u); ratmp = similar(u, uEltypeNoUnits)
  rk4cache = RK4Cache(u,uprev,rk1,rk2,rk3,rk4,rk,rtmp,ratmp)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(dt),5)
  c = fill(zero(t), 6, 6)
  g = fill(zero(t), 6)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 5)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 6)
  for i in 1:5
    ϕ_n[i] = zero(rate_prototype)
    ϕstar_nm1[i] = zero(rate_prototype)
    ϕstar_n[i] = zero(rate_prototype)
  end
  for i in 1:6
    ϕ_np1[i] = zero(rate_prototype)
  end
  β = fill(zero(t),5)
  order = 5
  atmp = similar(u,uEltypeNoUnits)
  tmp = zero(u)
  utilde = zero(u)
  VCABM5Cache(u,uprev,fsalfirst,rk4cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,atmp,tmp,utilde,1)
end

# VCABM

@cache mutable struct VCABMConstantCache{tArrayType,rArrayType,cArrayType,dtType,dtArrayType} <: OrdinaryDiffEqConstantCache
  ϕstar_nm1::rArrayType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::rArrayType
  ϕ_np1::rArrayType
  ϕstar_n::rArrayType
  β::tArrayType
  ξ::dtType
  ξ0::dtType
  order::Int
  max_order::Int
  step::Int
end

@cache mutable struct VCABMCache{uType,rateType,dtType,tArrayType,cArrayType,uNoUnitsType,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k4::rateType
  ϕstar_nm1::coefType
  dts::dtArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::coefType
  ϕ_np1::coefType
  ϕstar_n::coefType
  β::tArrayType
  order::Int
  max_order::Int
  atmp::uNoUnitsType
  tmp::uType
  ξ::dtType
  ξ0::dtType
  utilde::uType
  utildem1::uType
  utildem2::uType
  utildep1::uType
  atmpm1::uNoUnitsType
  atmpm2::uNoUnitsType
  atmpp1::uNoUnitsType
  step::Int
end

function alg_cache(alg::VCABM,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  dts = fill(zero(dt),13)
  c = fill(zero(t), 13, 13)
  g = fill(zero(t), 13)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 13)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 13)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 13)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 14)
  for i in 1:13
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(t),13)
  ξ = zero(dt)
  ξ0 = zero(dt)
  order = 1
  max_order = 12
  VCABMConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,ξ,ξ0,order,max_order,1)
end

function alg_cache(alg::VCABM,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(dt),13)
  c = fill(zero(t), 13, 13)
  g = fill(zero(t), 13)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 13)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 13)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 13)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 14)
  for i in 1:13
    ϕ_n[i] = zero(rate_prototype)
    ϕstar_nm1[i] = zero(rate_prototype)
    ϕstar_n[i] = zero(rate_prototype)
  end
  for i in 1:14
    ϕ_np1[i] = zero(rate_prototype)
  end
  β = fill(zero(t), 13)
  order = 1
  max_order = 12
  atmp = similar(u,uEltypeNoUnits)
  tmp = zero(u)
  ξ = zero(dt)
  ξ0 = zero(dt)
  utilde = zero(u)
  utildem2 = zero(u)
  utildem1 = zero(u)
  utildep1 = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  atmpm1 = similar(u,uEltypeNoUnits)
  atmpm2 = similar(u,uEltypeNoUnits)
  atmpp1 = similar(u,uEltypeNoUnits)
  VCABMCache(u,uprev,fsalfirst,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,max_order,atmp,tmp,ξ,ξ0,utilde,utildem1,utildem2,utildep1,atmpm1,atmpm2,atmpp1,1)
end

# IMEX Multistep methods

# CNAB2

@cache mutable struct CNAB2ConstantCache{rateType,N,uType,tType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  nlsolver::N
  uprev3::uType
  tprev2::tType
end

@cache mutable struct CNAB2Cache{uType,rateType,N,tType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  k1::rateType
  k2::rateType
  du₁::rateType
  nlsolver::N
  uprev3::uType
  tprev2::tType
end

function alg_cache(alg::CNAB2,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  γ, c = 1//2, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},γ,c,Val(false))

  k2 = rate_prototype
  uprev3 = u
  tprev2 = t

  CNAB2ConstantCache(k2,nlsolver,uprev3,tprev2)
end

function alg_cache(alg::CNAB2,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1//2, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  du₁ = zero(rate_prototype)
  uprev3 = zero(u)
  tprev2 = t

  CNAB2Cache(u,uprev,uprev2,fsalfirst,k1,k2,du₁,nlsolver,uprev3,tprev2)
end

# CNLF2

@cache mutable struct CNLF2ConstantCache{rateType,N,uType,tType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  nlsolver::N
  uprev2::uType
  uprev3::uType
  tprev2::tType
end

@cache mutable struct CNLF2Cache{uType,rateType,N,tType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  k1::rateType
  k2::rateType
  du₁::rateType
  nlsolver::N
  uprev3::uType
  tprev2::tType
end

function alg_cache(alg::CNLF2,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  γ, c = 1//1, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},γ,c,Val(false))

  k2 = rate_prototype
  uprev2 = u
  uprev3 = u
  tprev2 = t

  CNLF2ConstantCache(k2,nlsolver,uprev2,uprev3,tprev2)
end

function alg_cache(alg::CNLF2,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1//1, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  du₁ = zero(rate_prototype)
  uprev2 = zero(u)
  uprev3 = zero(u)
  tprev2 = t

  CNLF2Cache(u,uprev,uprev2,fsalfirst,k1,k2,du₁,nlsolver,uprev3,tprev2)
end
