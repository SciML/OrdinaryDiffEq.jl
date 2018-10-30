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
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  ralk2 = zero(rate_prototype)
  k  = zero(rate_prototype)
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
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  ralk2 = zero(rate_prototype)
  k  = zero(rate_prototype)
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
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  ralk2 = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = similar(u)
  t2 = zero(rate_prototype)
  t3 = zero(rate_prototype)
  t4 = zero(rate_prototype)
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
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  ralk2 = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = similar(u)
  t2 = zero(rate_prototype)
  t3 = zero(rate_prototype)
  t4 = zero(rate_prototype)
  t5 = zero(rate_prototype)
  t6 = zero(rate_prototype)
  t7 = zero(rate_prototype)
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
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = similar(u)
  t2 = zero(rate_prototype)
  t3 = zero(rate_prototype)
  t4 = zero(rate_prototype)
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
  fsalfirst = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = similar(u)
  t2 = zero(rate_prototype)
  t3 = zero(rate_prototype)
  t4 = zero(rate_prototype)
  t5 = zero(rate_prototype)
  t6 = zero(rate_prototype)
  t7 = zero(rate_prototype)
  t8 = zero(rate_prototype)
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
  dts = fill(zero(typeof(dt)),3)
  c = fill(zero(typeof(t)), 3, 3)
  g = fill(zero(typeof(t)), 3)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 3)
  for i in 1:3
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(typeof(t)),3)
  order = 3
  tab = BS3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  VCAB3ConstantCache(dts,c,g,ϕ_n,ϕstar_nm1,ϕstar_n,β,order,tab,1)
end

function alg_cache(alg::VCAB3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = BS3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  bk1 = zero(rate_prototype)
  bk2 = zero(rate_prototype)
  bk3 = zero(rate_prototype)
  bk4 = zero(rate_prototype)
  butilde = similar(u)
  batmp = similar(u,uEltypeNoUnits)
  btmp = similar(u)
  bs3cache = BS3Cache(u,uprev,bk1,bk2,bk3,bk4,butilde,btmp,batmp,tab)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(typeof(dt)),3)
  c = fill(zero(typeof(t)), 3, 3)
  g = fill(zero(typeof(t)), 3)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 3)
  for i in 1:3
    ϕ_n[i] = zero(rate_prototype)
    ϕstar_nm1[i] = zero(rate_prototype)
    ϕstar_n[i] = zero(rate_prototype)
  end
  β = fill(zero(typeof(t)),3)
  order = 3
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  utilde = similar(u)
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
  dts = fill(zero(typeof(dt)),4)
  c = fill(zero(typeof(t)), 4, 4)
  g = fill(zero(typeof(t)), 4)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 4)
  for i in 1:4
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(typeof(t)),4)
  order = 4
  rk4constcache = RK4ConstantCache()
  VCAB4ConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,rk4constcache,1)
end

function alg_cache(alg::VCAB4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  rk1 = zero(rate_prototype)
  rk2 = zero(rate_prototype)
  rk3 = zero(rate_prototype)
  rk4 = zero(rate_prototype)
  rk  = zero(rate_prototype)
  rtmp = similar(u); ratmp = similar(u, uEltypeNoUnits)
  rk4cache = RK4Cache(u,uprev,rk1,rk2,rk3,rk4,rk,rtmp,ratmp)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(typeof(dt)),4)
  c = fill(zero(typeof(t)), 4, 4)
  g = fill(zero(typeof(t)), 4)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 4)
  for i in 1:4
    ϕ_n[i] = zero(rate_prototype)
    ϕstar_nm1[i] = zero(rate_prototype)
    ϕstar_n[i] = zero(rate_prototype)
  end
  β = fill(zero(typeof(t)),4)
  order = 4
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  utilde = similar(u)
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
  dts = fill(zero(typeof(dt)),5)
  c = fill(zero(typeof(t)), 5, 5)
  g = fill(zero(typeof(t)), 5)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 5)
  for i in 1:5
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(typeof(t)),5)
  order = 5
  rk4constcache = RK4ConstantCache()
  VCAB5ConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕstar_n,β,order,rk4constcache,1)
end

function alg_cache(alg::VCAB5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  rk1 = zero(rate_prototype)
  rk2 = zero(rate_prototype)
  rk3 = zero(rate_prototype)
  rk4 = zero(rate_prototype)
  rk  = zero(rate_prototype)
  rtmp = similar(u); ratmp = similar(u, uEltypeNoUnits)
  rk4cache = RK4Cache(u,uprev,rk1,rk2,rk3,rk4,rk,rtmp,ratmp)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(typeof(dt)),5)
  c = fill(zero(typeof(t)), 5, 5)
  g = fill(zero(typeof(t)), 5)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 5)
  for i in 1:5
    ϕ_n[i] = zero(rate_prototype)
    ϕstar_nm1[i] = zero(rate_prototype)
    ϕstar_n[i] = zero(rate_prototype)
  end
  β = fill(zero(typeof(t)),5)
  order = 5
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  utilde = similar(u)
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
  dts = fill(zero(typeof(dt)),3)
  c = fill(zero(typeof(t)), 4, 4)
  g = fill(zero(typeof(t)), 4)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 3)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 3)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 4)
  for i in 1:3
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(typeof(t)),3)
  order = 3
  tab = BS3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  VCABM3ConstantCache(dts,c,g,ϕ_n,ϕ_np1,ϕstar_nm1,ϕstar_n,β,order,tab,1)
end

function alg_cache(alg::VCABM3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = BS3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  bk1 = zero(rate_prototype)
  bk2 = zero(rate_prototype)
  bk3 = zero(rate_prototype)
  bk4 = zero(rate_prototype)
  butilde = similar(u)
  batmp = similar(u,uEltypeNoUnits)
  btmp = similar(u)
  bs3cache = BS3Cache(u,uprev,bk1,bk2,bk3,bk4,butilde,btmp,batmp,tab)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(typeof(dt)),3)
  c = fill(zero(typeof(t)), 4, 4)
  g = fill(zero(typeof(t)), 4)
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
  β = fill(zero(typeof(t)),3)
  order = 3
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  utilde = similar(u)
  VCABM3Cache(u,uprev,fsalfirst,bs3cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,atmp,tmp,utilde,tab,1)
end

# VCABM4

mutable struct VCABM4ConstantCache{rk4constcache,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
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

mutable struct VCABM4Cache{uType,rateType,uArrayType,rk4cacheType,tArrayType,cArrayType,uEltypeNoUnits,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
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
  atmp::uEltypeNoUnits
  tmp::uType
  utilde::uArrayType
  step::Int
end

u_cache(c::VCABM4Cache) = ()
du_cache(c::VCABM4Cache) = ()

function alg_cache(alg::VCABM4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  dts = fill(zero(typeof(dt)),4)
  c = fill(zero(typeof(t)), 5, 5)
  g = fill(zero(typeof(t)), 5)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 4)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 4)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 5)
  for i in 1:4
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(typeof(t)),4)
  order = 4
  rk4constcache = RK4ConstantCache()
  VCABM4ConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,rk4constcache,1)
end

function alg_cache(alg::VCABM4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  rk1 = zero(rate_prototype)
  rk2 = zero(rate_prototype)
  rk3 = zero(rate_prototype)
  rk4 = zero(rate_prototype)
  rk  = zero(rate_prototype)
  rtmp = similar(u); ratmp = similar(u, uEltypeNoUnits)
  rk4cache = RK4Cache(u,uprev,rk1,rk2,rk3,rk4,rk,rtmp,ratmp)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(typeof(dt)),4)
  c = fill(zero(typeof(t)), 5, 5)
  g = fill(zero(typeof(t)), 5)
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
  β = fill(zero(typeof(t)),4)
  order = 4
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  utilde = similar(u)
  VCABM4Cache(u,uprev,fsalfirst,rk4cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,atmp,tmp,utilde,1)
end

# VCABM5

mutable struct VCABM5ConstantCache{rk4constcache,tArrayType,rArrayType,cArrayType,dtArrayType} <: OrdinaryDiffEqConstantCache
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

mutable struct VCABM5Cache{uType,rateType,uArrayType,rk4cacheType,tArrayType,cArrayType,uEltypeNoUnits,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
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
  atmp::uEltypeNoUnits
  tmp::uType
  utilde::uArrayType
  step::Int
end

u_cache(c::VCABM5Cache) = ()
du_cache(c::VCABM5Cache) = ()

function alg_cache(alg::VCABM5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  dts = fill(zero(typeof(dt)),5)
  c = fill(zero(typeof(t)), 6, 6)
  g = fill(zero(typeof(t)), 6)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 5)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 5)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 6)
  for i in 1:5
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(typeof(t)),5)
  order = 5
  rk4constcache = RK4ConstantCache()
  VCABM5ConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,rk4constcache,1)
end

function alg_cache(alg::VCABM5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  rk1 = zero(rate_prototype)
  rk2 = zero(rate_prototype)
  rk3 = zero(rate_prototype)
  rk4 = zero(rate_prototype)
  rk  = zero(rate_prototype)
  rtmp = similar(u); ratmp = similar(u, uEltypeNoUnits)
  rk4cache = RK4Cache(u,uprev,rk1,rk2,rk3,rk4,rk,rtmp,ratmp)
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(typeof(dt)),5)
  c = fill(zero(typeof(t)), 6, 6)
  g = fill(zero(typeof(t)), 6)
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
  β = fill(zero(typeof(t)),5)
  order = 5
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  utilde = similar(u)
  VCABM5Cache(u,uprev,fsalfirst,rk4cache,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,atmp,tmp,utilde,1)
end

# VCABM

mutable struct VCABMConstantCache{tArrayType,rArrayType,cArrayType,dtType,dtArrayType} <: OrdinaryDiffEqConstantCache
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

mutable struct VCABMCache{uType,rateType,uArrayType,dtType,tArrayType,cArrayType,uEltypeNoUnits,coefType,dtArrayType} <: OrdinaryDiffEqMutableCache
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
  atmp::uEltypeNoUnits
  tmp::uType
  ξ::dtType
  ξ0::dtType
  utilde::uArrayType
  utildem1::uArrayType
  utildem2::uArrayType
  utildep1::uArrayType
  atmpm1::uEltypeNoUnits
  atmpm2::uEltypeNoUnits
  atmpp1::uEltypeNoUnits
  step::Int
end

u_cache(c::VCABMCache) = ()
du_cache(c::VCABMCache) = ()

function alg_cache(alg::VCABM,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  dts = fill(zero(typeof(dt)),13)
  c = fill(zero(typeof(t)), 13, 13)
  g = fill(zero(typeof(t)), 13)
  ϕ_n = Vector{typeof(rate_prototype)}(undef, 13)
  ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 13)
  ϕstar_n = Vector{typeof(rate_prototype)}(undef, 13)
  ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 14)
  for i in 1:13
    ϕ_n[i] = copy(rate_prototype)
    ϕstar_nm1[i] = copy(rate_prototype)
    ϕstar_n[i] = copy(rate_prototype)
  end
  β = fill(zero(typeof(t)),13)
  ξ = zero(dt)
  ξ0 = zero(dt)
  order = 1
  max_order = 12
  VCABMConstantCache(ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,ξ,ξ0,order,max_order,1)
end

function alg_cache(alg::VCABM,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  fsalfirst = zero(rate_prototype)
  k4 = zero(rate_prototype)
  dts = fill(zero(typeof(dt)),13)
  c = fill(zero(typeof(t)), 13, 13)
  g = fill(zero(typeof(t)), 13)
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
  β = fill(zero(typeof(t)), 13)
  order = 1
  max_order = 12
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  ξ = zero(dt)
  ξ0 = zero(dt)
  utilde = similar(u)
  utildem2 = similar(u)
  utildem1 = similar(u)
  utildep1 = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  atmpm1 = similar(u,uEltypeNoUnits)
  atmpm2 = similar(u,uEltypeNoUnits)
  atmpp1 = similar(u,uEltypeNoUnits)
  VCABMCache(u,uprev,fsalfirst,k4,ϕstar_nm1,dts,c,g,ϕ_n,ϕ_np1,ϕstar_n,β,order,max_order,atmp,tmp,ξ,ξ0,utilde,utildem1,utildem2,utildep1,atmpm1,atmpm2,atmpp1,1)
end

# IMEX Multistep methods

# CNAB2

mutable struct CNAB2ConstantCache{rateType,F,N,uType,tType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  uf::F
  nlsolve::N
  uprev3::uType
  tprev2::tType
end

mutable struct CNAB2Cache{uType,rateType,uNoUnitsType,J,W,UF,JC,N,tType,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  k::rateType
  k1::rateType
  k2::rateType
  du₁::rateType
  du1::rateType
  z::uType
  dz::uType
  b::uType
  tmp::uType
  atmp::uNoUnitsType
  J::J
  W::W
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
  uprev3::uType
  tprev2::tType
end

u_cache(c::CNAB2Cache)    = ()
du_cache(c::CNAB2Cache)   = ()

function alg_cache(alg::CNAB2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  k2 = rate_prototype
  uprev3 = u
  tprev2 = t

  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1//2,1,ηold,z₊,dz,tmp,b,k))
  CNAB2ConstantCache(k2,uf,nlsolve,uprev3,tprev2)
end

function alg_cache(alg::CNAB2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  atmp = similar(u,uEltypeNoUnits)
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  du₁ = zero(rate_prototype)

  uprev3 = similar(u)
  tprev2 = t

  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1//2,1,ηold,z₊,dz,tmp,b,k))
  CNAB2Cache(u,uprev,uprev2,fsalfirst,k,k1,k2,du₁,du1,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,nlsolve,uprev3,tprev2)
end

# CNLF2

mutable struct CNLF2ConstantCache{rateType,F,N,uType,tType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  uf::F
  nlsolve::N
  uprev2::uType
  uprev3::uType
  tprev2::tType
end

mutable struct CNLF2Cache{uType,rateType,uNoUnitsType,J,W,UF,JC,N,tType,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  k::rateType
  k1::rateType
  k2::rateType
  du₁::rateType
  du1::rateType
  z::uType
  dz::uType
  b::uType
  tmp::uType
  atmp::uNoUnitsType
  J::J
  W::W
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
  uprev3::uType
  tprev2::tType
end

u_cache(c::CNLF2Cache)    = ()
du_cache(c::CNLF2Cache)   = ()

function alg_cache(alg::CNLF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  k2 = rate_prototype
  uprev2 = u
  uprev3 = u
  tprev2 = t

  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1//1,1,ηold,z₊,dz,tmp,b,k))
  CNLF2ConstantCache(k2,uf,nlsolve,uprev2,uprev3,tprev2)
end

function alg_cache(alg::CNLF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  atmp = similar(u,uEltypeNoUnits)
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  du₁ = zero(rate_prototype)

  uprev2 = similar(u)
  uprev3 = similar(u)
  tprev2 = t

  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1//1,1,ηold,z₊,dz,tmp,b,k))
  CNLF2Cache(u,uprev,uprev2,fsalfirst,k,k1,k2,du₁,du1,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,nlsolve,uprev3,tprev2)
end
