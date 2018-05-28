mutable struct AN5ConstantCache{zType,lType,dtType,uType,tsit5Type} <: OrdinaryDiffEqConstantCache
  # `z` is the Nordsieck vector
  z::zType
  # `l` is used for the corrector iteration
  l::MVector{6,lType}
  # `m` is a tmp vector that is used for calculating `l`
  m::MVector{6,lType}
  # `c_LTE` is used for the error estimation for the current order
  c_LTE::lType
  # `tau` stores `dt`s
  tau::MVector{6, dtType}
  # `Δ` is the difference between the predictor `uₙ₀` and `uₙ`
  Δ::uType
  # `Tsit5` for the first step
  tsit5tab::tsit5Type
  step::Int
end

function AN5ConstantCache(u, uprev, rate_prototype, uBottomEltypeNoUnits, tTypeNoUnits, dt)
  N = 5
  z = [zero(rate_prototype) for i in 1:N+1]
  Δ = u
  l = zeros(MVector{N+1,tTypeNoUnits}); m = zeros(l)
  c_LTE = zero(tTypeNoUnits)
  tau = zeros(MVector{N+1,typeof(dt)})
  tsit5tab = Tsit5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  AN5ConstantCache(z,l,m,c_LTE,tau,Δ,tsit5tab,1)
end

function alg_cache(alg::AN5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  AN5ConstantCache(u, uprev, rate_prototype, uBottomEltypeNoUnits, tTypeNoUnits, dt)
end

mutable struct AN5Cache{uType,rateType,uArrayType,uEltypeNoUnits,constType,tsit5Type} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  utilde::uArrayType
  tmp::uType
  ratetmp::rateType
  atmp::uEltypeNoUnits
  const_cache::constType
  tsit5cache::tsit5Type
end

u_cache(c::AN5Cache) = ()
du_cache(c::AN5Cache) = (c.fsalfirst,)

function alg_cache(alg::AN5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  fsalfirst = zeros(rate_prototype)
  const_cache = AN5ConstantCache(u, uprev, rate_prototype, uBottomEltypeNoUnits, tTypeNoUnits, dt)
  #################################################
  # Tsit5
  tab = const_cache.tsit5tab
  # Cannot alias pointers, since we have to use `k`s to start the Nordsieck vector
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype); k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype); k7 = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u)); tmp = similar(u)
  tsit5cache = Tsit5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp,tab)
  #################################################
  for i in 1:6
    const_cache.z[i] = zeros(rate_prototype)
  end
  ratetmp = k6
  const_cache.Δ = k7
  AN5Cache(u,uprev,fsalfirst,utilde,tmp,ratetmp,atmp,const_cache,tsit5cache)
end

mutable struct JVODEConstantCache{zType,lType,dtType,uType,tsit5Type} <: OrdinaryDiffEqConstantCache
  # `z` is the Nordsieck vector
  z::zType
  # `l` is used for the corrector iteration
  l::MVector{13,lType}
  # `m` is a tmp vector that is used for calculating `l`
  m::MVector{13,lType}
  # `c_LTE` is used for the error estimation for the current order + 1
  c_LTE₊₁::lType
  # `c_LTE` is used for the error estimation for the current order
  c_LTE::lType
  # `c_LTE` is used for the error estimation for the current order - 1
  c_LTE₋₁::lType
  # `tau` stores `dt`s
  tau::MVector{13, dtType}
  # `Δ` is the difference between the predictor `uₙ₀` and `uₙ`
  Δ::uType
  # `Tsit5` for the first step
  tsit5tab::tsit5Type
  step::Int
  nextorder::Int
  # number of steps to take before considering to change order
  n_wait::Int
end

function JVODEConstantCache(u, uprev, rate_prototype, uBottomEltypeNoUnits, tTypeNoUnits, dt)
  N = 12
  z = [rate_prototype for i in 1:N+1]
  Δ = u
  l = zeros(MVector{N+1,tTypeNoUnits}); m = zeros(l)
  c_LTE₊₁ = zero(tTypeNoUnits)
  c_LTE = zero(tTypeNoUnits)
  c_LTE₋₁ = zero(tTypeNoUnits)
  tau = zeros(MVector{N+1,typeof(dt)})
  tsit5tab = Tsit5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  JVODEConstantCache(z,l,m,c_LTE₊₁,c_LTE,c_LTE₋₁,tau,Δ,tsit5tab,4,4,4)
end

function alg_cache(alg::JVODE,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  JVODEConstantCache(u, uprev, rate_prototype, uBottomEltypeNoUnits, tTypeNoUnits, dt)
end

mutable struct JVODECache{uType,rateType,uArrayType,uEltypeNoUnits,constType,tsit5Type} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  utilde::uArrayType
  tmp::uType
  ratetmp::rateType
  atmp::uEltypeNoUnits
  const_cache::constType
  tsit5cache::tsit5Type
end

u_cache(c::JVODECache) = ()
du_cache(c::JVODECache) = (c.fsalfirst,)

function alg_cache(alg::JVODE,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  fsalfirst = zeros(rate_prototype)
  const_cache = JVODEConstantCache(u, uprev, rate_prototype, uBottomEltypeNoUnits, tTypeNoUnits, dt)
  #################################################
  # Tsit5
  tab = const_cache.tsit5tab
  # Cannot alias pointers, since we have to use `k`s to start the Nordsieck vector
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype); k7 = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u)); tmp = similar(u)
  tsit5cache = Tsit5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp,tab)
  #################################################
  # Nordsieck Vector
  z = const_cache.z
  # One-shot start
  z[1] = zeros(rate_prototype); z[2] = zeros(rate_prototype); z[3] = zeros(rate_prototype);
  z[4] = zeros(rate_prototype); z[5] = zeros(rate_prototype);
  # Order increase
  z[5] = k1; z[6] = k2; z[7] = k3; z[8] = k4; z[9] = k5; z[10] = k6; z[11] = k7;
  z[12] = zeros(rate_prototype); z[13] = zeros(rate_prototype);
  ratetmp = zeros(rate_prototype)
  #################################################
  const_cache.Δ = utilde
  JVODECache(u,uprev,fsalfirst,utilde,tmp,ratetmp,atmp,const_cache,tsit5cache)
end
