mutable struct AN5ConstantCache{zType,lType,dtType,uType,vern9Type} <: OrdinaryDiffEqConstantCache
  # `z` is the Nordsieck vector
  z::zType
  # `l` is used for the corrector iteration
  l::MVector{6,lType}
  # `m` is a tmp vector that is used for calculating `l`
  m::MVector{6,lType}
  # `tq` is used for the error estimation for the current order
  tq::lType
  # `tau` stores `dt`s
  tau::MVector{6, dtType}
  # `Δ` is the difference between the predictor `uₙ₀` and `uₙ`
  Δ::uType
  # `Tsit5` for the first step
  vern9tab::vern9Type
  step::Int
end

function AN5ConstantCache(u, uprev, rate_prototype, uBottomEltypeNoUnits, tTypeNoUnits, dt, mut)
  #siz = size(rate_prototype)
  #typ = Base.promote_op(*, eltype(rate_prototype), tTypeNoUnits)
  if mut
    # We don't need to swap pointers for the mutable cache
    z = SVector(zeros(u), [zeros(rate_prototype) for i in 1:5]...)
    Δ = zeros(u)
  else
    z = [u, [copy(rate_prototype) for i in 1:5]...]
    Δ = u
  end
  l = zeros(MVector{6,tTypeNoUnits}); m = zeros(l)
  tq = zero(tTypeNoUnits)
  tau = zeros(MVector{6,typeof(dt)})
  vern9tab = Vern9ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  AN5ConstantCache(z,l,m,tq,tau,Δ,vern9tab,1)
end

function alg_cache(alg::AN5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  AN5ConstantCache(u, uprev, rate_prototype, uBottomEltypeNoUnits, tTypeNoUnits, dt, false)
end

mutable struct AN5Cache{uType,rateType,uArrayType,uEltypeNoUnits,constType,vern9Type} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  utilde::uArrayType
  tmp::uType
  ratetmp::rateType
  atmp::uEltypeNoUnits
  const_cache::constType
  vern9cache::vern9Type
end

u_cache(c::AN5Cache) = ()
du_cache(c::AN5Cache) = (c.fsalfirst,)

function alg_cache(alg::AN5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  fsalfirst = zeros(rate_prototype)
  const_cache = AN5ConstantCache(u, uprev, rate_prototype, uBottomEltypeNoUnits, tTypeNoUnits, dt, true)
  #################################################
  # Tsit5
  tab = const_cache.vern9tab
  #k1 = const_cache.z[2]; k2 = const_cache.z[3]
  #k3 = const_cache.z[4]; k4 = const_cache.z[5]
  #k5 = const_cache.z[6]
  # Cannot alias pointers, since we have to use `k`s to start the Nordsieck vector
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype); k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype); k7 = zeros(rate_prototype)
  k8 = zeros(rate_prototype); k9 = zeros(rate_prototype)
  k10 = zeros(rate_prototype); k11 = zeros(rate_prototype)
  k12 = zeros(rate_prototype); k13 = zeros(rate_prototype)
  k14 = zeros(rate_prototype); k15 = zeros(rate_prototype)
  k16 = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u)); tmp = similar(u)
  vern9cache = Vern9Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,utilde,tmp,atmp,tab)
  #################################################
  ratetmp = k6
  AN5Cache(u,uprev,fsalfirst,utilde,tmp,ratetmp,atmp,const_cache,vern9cache)
end
