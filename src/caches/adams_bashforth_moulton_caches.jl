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

mutable struct VCAB3ConstantCache{rateType,TabType,tArrayType,rArrayType,cArrayType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  k3::rateType
  ϕstar_nm1::rArrayType
  grid_points::tArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::rArrayType
  ϕstar_n::rArrayType
  β::tArrayType
  k::Int
  order::Int
  tab::TabType
end

mutable struct VCAB3Cache{uType,rateType,TabType,uArrayType,bs3Type,tArrayType,cArrayType,uEltypeNoUnits,coefType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  bs3cache::bs3Type
  k2::rateType
  k3::rateType
  k4::rateType
  ϕstar_nm1::coefType
  grid_points::tArrayType
  c::cArrayType
  g::tArrayType
  ϕ_n::coefType
  ϕstar_n::coefType
  β::tArrayType
  k::Int
  order::Int
  atmp::uEltypeNoUnits
  tmp::uType
  utilde::uArrayType
  tab::TabType
end

u_cache(c::VCAB3Cache) = ()
du_cache(c::VCAB3Cache) = ()

function alg_cache(alg::VCAB3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  k2 = rate_prototype
  k3 = rate_prototype
  ϕstar_nm1 = zeros(typeof(rate_prototype),3)
  grid_points = zeros(typeof(t),3)
  c = zeros(typeof(t), 3, 3)
  g = zeros(typeof(t), 3)
  ϕ_n = zeros(typeof(rate_prototype),3)
  ϕstar_n = zeros(typeof(rate_prototype),3)
  β = zeros(typeof(t),3)
  k = 1
  order = 3
  tab = BS3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  VCAB3ConstantCache(k2,k3,ϕstar_nm1,grid_points,c,g,ϕ_n,ϕstar_n,β,k,order,tab)
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
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  ϕstar_nm1 = Array{typeof(rate_prototype)}(1,3)
  grid_points = zeros(typeof(t),3)
  c = zeros(typeof(t),3,3)
  g = zeros(typeof(t),3)
  ϕ_n = Array{typeof(rate_prototype)}(1,3)
  ϕstar_n = Array{typeof(rate_prototype)}(1,3)
  β = zeros(typeof(t),3)
  k = 1
  order = 3
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  utilde = similar(u,indices(u))
  VCAB3Cache(u,uprev,fsalfirst,bs3cache,k2,k3,k4,ϕstar_nm1,grid_points,c,g,ϕ_n,ϕstar_n,β,k,order,atmp,tmp,utilde,tab)
end
