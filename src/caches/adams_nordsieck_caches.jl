mutable struct AN5ConstantCache{zType,lType,dtType,uType} <: OrdinaryDiffEqConstantCache
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
  step::Int
end

function AN5ConstantCache(u, rate_prototype, tTypeNoUnits, dt, mut)
  #siz = size(rate_prototype)
  #typ = Base.promote_op(*, eltype(rate_prototype), tTypeNoUnits)
  if mut
    # We don't need to swap pointers for the mutable cache
    z = SVector(u, [zeros(rate_prototype) for i in 1:5]...)
    Δ = zeros(u)
  else
    z = [u, [rate_prototype for i in 1:5]...]
    Δ = u
  end
  l = zeros(tTypeNoUnits, MVector{6}); m = zeros(l)
  tq = zero(tTypeNoUnits)
  tau = zeros(dt, MVector{6})
  AN5ConstantCache(z,l,m,tq,tau,Δ,1)
end

function alg_cache(alg::AN5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  AN5ConstantCache(u, rate_prototype, tTypeNoUnits, dt, false)
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
  tab = AN5ConstantCache(u, rate_prototype, tTypeNoUnits, dt, true)
  AN5Cache(u,uprev,fsalfirst,tmp,tab)
end
