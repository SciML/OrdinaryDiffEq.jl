mutable struct AB3Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
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

u_cache(c::AB3Cache) = ()
du_cache(c::AB3Cache) = (c.fsalfirst,c.k2,c.k3,c.ralk2,c.k)

mutable struct AB3ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  step::Int
end

function alg_cache(alg::AB3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  fsalfirst = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  ralk2 = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  AB3Cache(u,uprev,fsalfirst,k2,k3,ralk2,k,tmp,1)
end

function alg_cache(alg::AB3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  k2 = rate_prototype
  k3 = rate_prototype
  AB3ConstantCache(k2,k3,1)
end

mutable struct ABM32Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
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

u_cache(c::ABM32Cache) = ()
du_cache(c::ABM32Cache) = (c.fsalfirst,c.k2,c.k3,c.ralk2,c.k)

mutable struct ABM32ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  step::Int
end

function alg_cache(alg::ABM32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  fsalfirst = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  ralk2 = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  ABM32Cache(u,uprev,fsalfirst,k2,k3,ralk2,k,tmp,1)
end

function alg_cache(alg::ABM32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  k2 = rate_prototype
  k3 = rate_prototype
  ABM32ConstantCache(k2,k3,1)
end

mutable struct AB4Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
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

u_cache(c::AB4Cache) = ()
du_cache(c::AB4Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.ralk2,c.k)

mutable struct AB4ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  k4::rateType
  step::Int
end

function alg_cache(alg::AB4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  fsalfirst = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  ralk2 = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  t2 = zeros(rate_prototype)
  t3 = zeros(rate_prototype)
  t4 = zeros(rate_prototype)
  AB4Cache(u,uprev,fsalfirst,k2,k3,k4,ralk2,k,tmp,t2,t3,t4,1)
end

function alg_cache(alg::AB4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  k2 = rate_prototype
  k3 = rate_prototype
  k4 = rate_prototype
  AB4ConstantCache(k2,k3,k4,1)
end

mutable struct ABM43Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
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

u_cache(c::ABM43Cache) = ()
du_cache(c::ABM43Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.ralk2,c.k)

mutable struct ABM43ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  k4::rateType
  step::Int
end

function alg_cache(alg::ABM43,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  fsalfirst = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  ralk2 = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  t2 = zeros(rate_prototype)
  t3 = zeros(rate_prototype)
  t4 = zeros(rate_prototype)
  t5 = zeros(rate_prototype)
  t6 = zeros(rate_prototype)
  t7 = zeros(rate_prototype)
  ABM43Cache(u,uprev,fsalfirst,k2,k3,k4,ralk2,k,tmp,t2,t3,t4,t5,t6,t7,1)
end

function alg_cache(alg::ABM43,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  k2 = rate_prototype
  k3 = rate_prototype
  k4 = rate_prototype
  ABM43ConstantCache(k2,k3,k4,1)
end

mutable struct AB5Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
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

u_cache(c::AB5Cache) = ()
du_cache(c::AB5Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k)

mutable struct AB5ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  step::Int
end

function alg_cache(alg::AB5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  fsalfirst = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  t2 = zeros(rate_prototype)
  t3 = zeros(rate_prototype)
  t4 = zeros(rate_prototype)
  AB5Cache(u,uprev,fsalfirst,k2,k3,k4,k5,k,tmp,t2,t3,t4,1)
end

function alg_cache(alg::AB5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  k2 = rate_prototype
  k3 = rate_prototype
  k4 = rate_prototype
  k5 = rate_prototype
  AB5ConstantCache(k2,k3,k4,k5,1)
end

mutable struct ABM54Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
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

u_cache(c::ABM54Cache) = ()
du_cache(c::ABM54Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k)

mutable struct ABM54ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  step::Int
end

function alg_cache(alg::ABM54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  fsalfirst = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  t2 = zeros(rate_prototype)
  t3 = zeros(rate_prototype)
  t4 = zeros(rate_prototype)
  t5 = zeros(rate_prototype)
  t6 = zeros(rate_prototype)
  t7 = zeros(rate_prototype)
  t8 = zeros(rate_prototype)
  ABM54Cache(u,uprev,fsalfirst,k2,k3,k4,k5,k,tmp,t2,t3,t4,t5,t6,t7,t8,1)
end

function alg_cache(alg::ABM54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  k2 = rate_prototype
  k3 = rate_prototype
  k4 = rate_prototype
  k5 = rate_prototype
  ABM54ConstantCache(k2,k3,k4,k5,1)
end

mutable struct VCAB3ConstantCache{TabType,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
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

mutable struct VCAB3Cache{uType,rateType,TabType,uArrayType,bs3Type,tArrayType,cArrayType,uEltypeNoUnits,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
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
  atmp::uEltypeNoUnits
  tmp::uType
  utilde::uArrayType
  tab::TabType
  step::Int
end

u_cache(c::VCAB3Cache) = ()
du_cache(c::VCAB3Cache) = ()

function alg_cache(alg::VCAB3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  dts = zeros(typeof(dt),3)
  c = zeros(typeof(t), 3, 3)
  g = zeros(typeof(t), 3)
  ϕ_n = Vector{typeof(rate_prototype)}(3)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(3)
  ϕstar_n = Vector{typeof(rate_prototype)}(3)
  for i in 1:3
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = zeros(typeof(t),3)
  order = 3
  tab = BS3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  VCAB3ConstantCache(dts,c,g,ϕ_n,ϕstar_nm1,ϕstar_n,β,order,tab,1)
end

function alg_cache(alg::VCAB3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = BS3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  bk1 = zeros(rate_prototype)
  bk2 = zeros(rate_prototype)
  bk3 = zeros(rate_prototype)
  bk4 = zeros(rate_prototype)
  butilde = similar(u,indices(u))
  batmp = similar(u,uEltypeNoUnits)
  btmp = similar(u)
  bs3cache = BS3Cache(u,uprev,bk1,bk2,bk3,bk4,butilde,btmp,batmp,tab)
  fsalfirst = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  dts = zeros(typeof(dt),3)
  c = zeros(typeof(t),3,3)
  g = zeros(typeof(t),3)
  ϕ_n = Vector{typeof(rate_prototype)}(3)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(3)
  ϕstar_n = Vector{typeof(rate_prototype)}(3)
  for i in 1:3
    ϕ_n[i] = zeros(rate_prototype)
    ϕstar_nm1[i] = zeros(rate_prototype)
    ϕstar_n[i] = zeros(rate_prototype)
  end
  β = zeros(typeof(t),3)
  order = 3
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  utilde = similar(u,indices(u))
  VCAB3Cache(u,uprev,fsalfirst,bs3cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,atmp,tmp,utilde,tab,1)
end

mutable struct VCAB4ConstantCache{rk4constcache,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
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

mutable struct VCAB4Cache{uType,rateType,uArrayType,rk4cacheType,tArrayType,cArrayType,uEltypeNoUnits,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
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
  atmp::uEltypeNoUnits
  tmp::uType
  utilde::uArrayType
  step::Int
end

u_cache(c::VCAB4Cache) = ()
du_cache(c::VCAB4Cache) = ()

function alg_cache(alg::VCAB4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  dts = zeros(typeof(dt),4)
  c = zeros(typeof(t), 4, 4)
  g = zeros(typeof(t), 4)
  ϕ_n = Vector{typeof(rate_prototype)}(4)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(4)
  ϕstar_n = Vector{typeof(rate_prototype)}(4)
  for i in 1:4
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = zeros(typeof(t),4)
  order = 4
  rk4constcache = RK4ConstantCache()
  VCAB4ConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,rk4constcache,1)
end

function alg_cache(alg::VCAB4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  rk1 = zeros(rate_prototype)
  rk2 = zeros(rate_prototype)
  rk3 = zeros(rate_prototype)
  rk4 = zeros(rate_prototype)
  rk  = zeros(rate_prototype)
  rtmp = similar(u); ratmp = similar(u, uEltypeNoUnits)
  rk4cache = RK4Cache(u,uprev,rk1,rk2,rk3,rk4,rk,rtmp,ratmp)
  fsalfirst = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  dts = zeros(typeof(dt),4)
  c = zeros(typeof(t),4,4)
  g = zeros(typeof(t),4)
  ϕ_n = Vector{typeof(rate_prototype)}(4)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(4)
  ϕstar_n = Vector{typeof(rate_prototype)}(4)
  for i in 1:4
    ϕ_n[i] = zeros(rate_prototype)
    ϕstar_nm1[i] = zeros(rate_prototype)
    ϕstar_n[i] = zeros(rate_prototype)
  end
  β = zeros(typeof(t),4)
  order = 4
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  utilde = similar(u,indices(u))
  VCAB4Cache(u,uprev,fsalfirst,rk4cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,atmp,tmp,utilde,1)
end

# VCAB5

mutable struct VCAB5ConstantCache{rk4constcache,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
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

mutable struct VCAB5Cache{uType,rateType,uArrayType,rk4cacheType,tArrayType,cArrayType,uEltypeNoUnits,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
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
  atmp::uEltypeNoUnits
  tmp::uType
  utilde::uArrayType
  step::Int
end

u_cache(c::VCAB5Cache) = ()
du_cache(c::VCAB5Cache) = ()

function alg_cache(alg::VCAB5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  dts = zeros(typeof(dt),5)
  c = zeros(typeof(t), 5, 5)
  g = zeros(typeof(t), 5)
  ϕ_n = Vector{typeof(rate_prototype)}(5)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(5)
  ϕstar_n = Vector{typeof(rate_prototype)}(5)
  for i in 1:5
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = zeros(typeof(t),5)
  order = 5
  rk4constcache = RK4ConstantCache()
  VCAB5ConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,rk4constcache,1)
end

function alg_cache(alg::VCAB5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  rk1 = zeros(rate_prototype)
  rk2 = zeros(rate_prototype)
  rk3 = zeros(rate_prototype)
  rk4 = zeros(rate_prototype)
  rk  = zeros(rate_prototype)
  rtmp = similar(u); ratmp = similar(u, uEltypeNoUnits)
  rk4cache = RK4Cache(u,uprev,rk1,rk2,rk3,rk4,rk,rtmp,ratmp)
  fsalfirst = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  dts = zeros(typeof(dt),5)
  c = zeros(typeof(t),5,5)
  g = zeros(typeof(t),5)
  ϕ_n = Vector{typeof(rate_prototype)}(5)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(5)
  ϕstar_n = Vector{typeof(rate_prototype)}(5)
  for i in 1:5
    ϕ_n[i] = zeros(rate_prototype)
    ϕstar_nm1[i] = zeros(rate_prototype)
    ϕstar_n[i] = zeros(rate_prototype)
  end
  β = zeros(typeof(t),5)
  order = 5
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  utilde = similar(u,indices(u))
  VCAB5Cache(u,uprev,fsalfirst,rk4cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,atmp,tmp,utilde,1)
end

# VCABM3

mutable struct VCABM3ConstantCache{TabType,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
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

mutable struct VCABM3Cache{uType,rateType,TabType,uArrayType,bs3Type,tArrayType,cArrayType,uEltypeNoUnits,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
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
  atmp::uEltypeNoUnits
  tmp::uType
  utilde::uArrayType
  tab::TabType
  step::Int
end

u_cache(c::VCABM3Cache) = ()
du_cache(c::VCABM3Cache) = ()

function alg_cache(alg::VCABM3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  dts = zeros(typeof(dt),3)
  c = zeros(typeof(t), 4, 4)
  g = zeros(typeof(t), 4)
  ϕ_n = Vector{typeof(rate_prototype)}(3)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(3)
  ϕstar_n = Vector{typeof(rate_prototype)}(3)
  ϕ_np1 = Vector{typeof(rate_prototype)}(4)
  for i in 1:3
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = zeros(typeof(t),3)
  order = 3
  tab = BS3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  VCABM3ConstantCache(dts,c,g,ϕ_n,ϕ_np1,ϕstar_nm1,ϕstar_n,β,order,tab,1)
end

function alg_cache(alg::VCABM3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = BS3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  bk1 = zeros(rate_prototype)
  bk2 = zeros(rate_prototype)
  bk3 = zeros(rate_prototype)
  bk4 = zeros(rate_prototype)
  butilde = similar(u,indices(u))
  batmp = similar(u,uEltypeNoUnits)
  btmp = similar(u)
  bs3cache = BS3Cache(u,uprev,bk1,bk2,bk3,bk4,butilde,btmp,batmp,tab)
  fsalfirst = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  dts = zeros(typeof(dt),3)
  c = zeros(typeof(t),4,4)
  g = zeros(typeof(t),4)
  ϕ_n = Vector{typeof(rate_prototype)}(3)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(3)
  ϕstar_n = Vector{typeof(rate_prototype)}(3)
  ϕ_np1 = Vector{typeof(rate_prototype)}(4)
  for i in 1:3
    ϕ_n[i] = zeros(rate_prototype)
    ϕstar_nm1[i] = zeros(rate_prototype)
    ϕstar_n[i] = zeros(rate_prototype)
  end
  β = zeros(typeof(t),3)
  order = 3
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  utilde = similar(u,indices(u))
  VCABM3Cache(u,uprev,fsalfirst,bs3cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,atmp,tmp,utilde,tab,1)
end
