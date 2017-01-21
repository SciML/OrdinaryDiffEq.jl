abstract OrdinaryDiffEqCache <: DECache
abstract OrdinaryDiffEqConstantCache <: OrdinaryDiffEqCache
abstract OrdinaryDiffEqMutableCache <: OrdinaryDiffEqCache
immutable ODEEmptyCache <: OrdinaryDiffEqConstantCache end
immutable ODEChunkCache{CS} <: OrdinaryDiffEqConstantCache end

type CompositeCache{T,F} <: OrdinaryDiffEqCache
  caches::T
  choice_function::F
  current::Int
end

function alg_cache{T}(alg::CompositeAlgorithm,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{T}})
  caches = map((x)->alg_cache(x,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,Val{T}),alg.algs)
  CompositeCache(caches,alg.choice_function,1)
end

alg_cache{F}(alg::OrdinaryDiffEqAlgorithm,prob,callback::F) = ODEEmptyCache()

immutable EulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  fsalfirst::rateType
end

function alg_cache(alg::Euler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  EulerCache(u,uprev,similar(rate_prototype,indices(rate_prototype)),similar(rate_prototype,indices(rate_prototype)))
end

immutable EulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Euler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = EulerConstantCache()

immutable MidpointCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  du::rateType
  utilde::uType
  fsalfirst::rateType
end

immutable MidpointConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::Midpoint,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  utilde = similar(u)
  k = similar(rate_prototype,indices(rate_prototype))
  du = similar(rate_prototype,indices(rate_prototype))
  fsalfirst = similar(rate_prototype,indices(rate_prototype))
  MidpointCache(u,uprev,k,du,utilde,fsalfirst)
end

alg_cache(alg::Midpoint,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = MidpointConstantCache()

immutable RK4Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k₂::rateType
  k₃::rateType
  k₄::rateType
  k::rateType
  tmp::uType
end

immutable RK4ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::RK4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  k₁ = similar(rate_prototype,indices(rate_prototype))
  k₂ = similar(rate_prototype,indices(rate_prototype))
  k₃ = similar(rate_prototype,indices(rate_prototype))
  k₄ = similar(rate_prototype,indices(rate_prototype))
  k  = similar(rate_prototype,indices(rate_prototype))
  tmp = similar(u)
  RK4Cache(u,uprev,k₁,k₂,k₃,k₄,k,tmp)
end

alg_cache(alg::RK4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = RK4ConstantCache()

immutable BS3Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

function alg_cache(alg::BS3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = BS3ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype))
  k2 = similar(rate_prototype,indices(rate_prototype))
  k3 = similar(rate_prototype,indices(rate_prototype))
  k4 = similar(rate_prototype,indices(rate_prototype))
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  BS3Cache(u,uprev,k1,k2,k3,k4,utilde,tmp,atmp,tab)
end

alg_cache(alg::BS3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = BS3ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable BS5Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  utilde::uType
  uhat::uType
  tmp::uType
  atmp::uEltypeNoUnits
  atmptilde::uEltypeNoUnits
  tab::TabType
end

function alg_cache(alg::BS5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = BS5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype))
  k2 = similar(rate_prototype,indices(rate_prototype))
  k3 = similar(rate_prototype,indices(rate_prototype))
  k4 = similar(rate_prototype,indices(rate_prototype))
  k5 = similar(rate_prototype,indices(rate_prototype))
  k6 = similar(rate_prototype,indices(rate_prototype))
  k7 = similar(rate_prototype,indices(rate_prototype))
  k8 = similar(rate_prototype,indices(rate_prototype))
  utilde = similar(u); uhat = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  atmptilde = similar(u,uEltypeNoUnits)
  BS5Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,utilde,uhat,tmp,atmp,atmptilde,tab)
end

alg_cache(alg::BS5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = BS5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Tsit5Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

function alg_cache(alg::Tsit5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = Tsit5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype))
  k2 = similar(rate_prototype,indices(rate_prototype))
  k3 = similar(rate_prototype,indices(rate_prototype))
  k4 = similar(rate_prototype,indices(rate_prototype))
  k5 = similar(rate_prototype,indices(rate_prototype))
  k6 = similar(rate_prototype,indices(rate_prototype))
  k7 = similar(rate_prototype,indices(rate_prototype))
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  atmptilde = similar(u,uEltypeNoUnits)
  Tsit5Cache(u,k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp,tab)
end

alg_cache(alg::Tsit5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = Tsit5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable DP5Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  dense_tmp3::rateType
  dense_tmp4::rateType
  update::rateType
  bspl::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

function alg_cache(alg::DP5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  k1 = similar(rate_prototype,indices(rate_prototype))
  k2 = similar(rate_prototype,indices(rate_prototype))
  k3 = similar(rate_prototype,indices(rate_prototype))
  k4 = similar(rate_prototype,indices(rate_prototype))
  k5 = similar(rate_prototype,indices(rate_prototype))
  k6 = similar(rate_prototype,indices(rate_prototype))
  k7 = similar(rate_prototype,indices(rate_prototype))
  dense_tmp3 = similar(rate_prototype,indices(rate_prototype))
  dense_tmp4 = similar(rate_prototype,indices(rate_prototype))
  update = similar(rate_prototype,indices(rate_prototype))
  bspl = similar(rate_prototype,indices(rate_prototype))
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  tab = DP5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  DP5Cache(u,k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp,tab)
end

alg_cache(alg::DP5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = DP5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable DP5ThreadedCache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  dense_tmp3::rateType
  dense_tmp4::rateType
  update::rateType
  bspl::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

function alg_cache(alg::DP5Threaded,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  k1 = similar(rate_prototype,indices(rate_prototype))
  k2 = similar(rate_prototype,indices(rate_prototype))
  k3 = similar(rate_prototype,indices(rate_prototype))
  k4 = similar(rate_prototype,indices(rate_prototype))
  k5 = similar(rate_prototype,indices(rate_prototype))
  k6 = similar(rate_prototype,indices(rate_prototype))
  k7 = similar(rate_prototype,indices(rate_prototype))
  dense_tmp3 = similar(rate_prototype,indices(rate_prototype))
  dense_tmp4 = similar(rate_prototype,indices(rate_prototype))
  update = similar(rate_prototype,indices(rate_prototype))
  bspl = similar(rate_prototype,indices(rate_prototype))
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  tab = DP5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  DP5ThreadedCache(u,k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp,tab)
end

immutable Vern6Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

function alg_cache(alg::Vern6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = Vern6ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype))
  k2 = similar(rate_prototype,indices(rate_prototype)); k3 = similar(rate_prototype,indices(rate_prototype)); k4 = similar(rate_prototype,indices(rate_prototype));
  k5 = similar(rate_prototype,indices(rate_prototype)); k6 = similar(rate_prototype,indices(rate_prototype)); k7 = similar(rate_prototype,indices(rate_prototype));
  k8 = similar(rate_prototype,indices(rate_prototype)); k9 = similar(rate_prototype,indices(rate_prototype));
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits);
  Vern6Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,utilde,tmp,atmp,tab)
end

alg_cache(alg::Vern6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = Vern6ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Vern7Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

function alg_cache(alg::Vern7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = Vern7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype)); k2 = similar(rate_prototype,indices(rate_prototype)); k3 = similar(rate_prototype,indices(rate_prototype)); k4 = similar(rate_prototype,indices(rate_prototype));
  k5 = similar(rate_prototype,indices(rate_prototype)); k6 = similar(rate_prototype,indices(rate_prototype)); k7 = similar(rate_prototype,indices(rate_prototype)); k8 = similar(rate_prototype,indices(rate_prototype));
  k9 = similar(rate_prototype,indices(rate_prototype)); k10 = similar(rate_prototype,indices(rate_prototype)); utilde = similar(u); update = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  Vern7Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = Vern7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))


immutable Vern8Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

function alg_cache(alg::Vern8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = Vern8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype)); k2 = similar(rate_prototype,indices(rate_prototype)); k3 = similar(rate_prototype,indices(rate_prototype));
  k4 = similar(rate_prototype,indices(rate_prototype));
  k5 = similar(rate_prototype,indices(rate_prototype)); k6 = similar(rate_prototype,indices(rate_prototype)); k7 = similar(rate_prototype,indices(rate_prototype));
  k8 = similar(rate_prototype,indices(rate_prototype)); tmp = similar(u)
  k9 = similar(rate_prototype,indices(rate_prototype)); k10 = similar(rate_prototype,indices(rate_prototype)); k11 = similar(rate_prototype,indices(rate_prototype));
  k12 = similar(rate_prototype,indices(rate_prototype)); k13 = similar(rate_prototype,indices(rate_prototype))
  utilde = similar(u); update = similar(u);
  atmp = similar(u,uEltypeNoUnits)
  Vern8Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = Vern8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Vern9Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

function alg_cache(alg::Vern9,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = Vern9ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype)); k2 = similar(rate_prototype,indices(rate_prototype));k3 = similar(rate_prototype,indices(rate_prototype));
  k4 = similar(rate_prototype,indices(rate_prototype));
  k5 = similar(rate_prototype,indices(rate_prototype)); k6 = similar(rate_prototype,indices(rate_prototype));k7 = similar(rate_prototype,indices(rate_prototype));
  k8 = similar(rate_prototype,indices(rate_prototype));
  k9 = similar(rate_prototype,indices(rate_prototype)); k10 = similar(rate_prototype,indices(rate_prototype)); k11 = similar(rate_prototype,indices(rate_prototype));
  k12 = similar(rate_prototype,indices(rate_prototype)); update = similar(u)
  k13 = similar(rate_prototype,indices(rate_prototype)); k14 = similar(rate_prototype,indices(rate_prototype)); k15 = similar(rate_prototype,indices(rate_prototype));
  k16 =similar(rate_prototype,indices(rate_prototype));
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits);
  Vern9Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern9,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = Vern9ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable TanYam7Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

function alg_cache(alg::TanYam7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = TanYam7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype)); k2 = similar(rate_prototype,indices(rate_prototype)) ; k3 = similar(rate_prototype,indices(rate_prototype)); k4 = similar(rate_prototype,indices(rate_prototype))
  k5 = similar(rate_prototype,indices(rate_prototype)); k6 = similar(rate_prototype,indices(rate_prototype)) ; k7 = similar(rate_prototype,indices(rate_prototype)); k8 = similar(rate_prototype,indices(rate_prototype))
  k9 = similar(rate_prototype,indices(rate_prototype)); k10= similar(rate_prototype,indices(rate_prototype)) ;
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype,indices(rate_prototype))
  TanYam7Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,tmp,atmp,k,tab)
end

alg_cache(alg::TanYam7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = TanYam7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))


immutable DP8Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  kupdate::rateType
  update::uType
  udiff::rateType
  bspl::rateType
  dense_tmp3::rateType
  dense_tmp4::rateType
  dense_tmp5::rateType
  dense_tmp6::rateType
  dense_tmp7::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  atmp2::uEltypeNoUnits
  tab::TabType
end

function alg_cache(alg::DP8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  k1 = similar(rate_prototype,indices(rate_prototype)); k2  = similar(rate_prototype,indices(rate_prototype)); k3  = similar(rate_prototype,indices(rate_prototype));  k4 = similar(rate_prototype,indices(rate_prototype))
  k5 = similar(rate_prototype,indices(rate_prototype)); k6  = similar(rate_prototype,indices(rate_prototype)); k7  = similar(rate_prototype,indices(rate_prototype));  k8 = similar(rate_prototype,indices(rate_prototype))
  k9 = similar(rate_prototype,indices(rate_prototype)); k10 = similar(rate_prototype,indices(rate_prototype)); k11 = similar(rate_prototype,indices(rate_prototype)); k12 = similar(rate_prototype,indices(rate_prototype))
  kupdate = similar(rate_prototype,indices(rate_prototype)); utilde = similar(u);
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); atmp2 = similar(u,uEltypeNoUnits); update = similar(u)
  k13 = similar(rate_prototype,indices(rate_prototype))
  k14 = similar(rate_prototype,indices(rate_prototype))
  k15 = similar(rate_prototype,indices(rate_prototype))
  k16 = similar(rate_prototype,indices(rate_prototype))
  udiff = similar(rate_prototype,indices(rate_prototype))
  bspl = similar(rate_prototype,indices(rate_prototype))
  # dense_tmp1 = udiff
  # dense_tmp2 = bspl
  dense_tmp3 = similar(rate_prototype,indices(rate_prototype))
  dense_tmp4 = similar(rate_prototype,indices(rate_prototype))
  dense_tmp5 = similar(rate_prototype,indices(rate_prototype))
  dense_tmp6 = similar(rate_prototype,indices(rate_prototype))
  dense_tmp7 = similar(rate_prototype,indices(rate_prototype))
  tab = DP8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  DP8Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,kupdate,
           update,udiff,bspl,dense_tmp3,dense_tmp4,dense_tmp5,dense_tmp6,dense_tmp7,
           utilde,tmp,atmp,atmp2,tab)
end

alg_cache(alg::DP8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = DP8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable TsitPap8Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

function alg_cache(alg::TsitPap8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = TsitPap8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype)); k2 = similar(rate_prototype,indices(rate_prototype)); k3 = similar(rate_prototype,indices(rate_prototype)); k4 = similar(rate_prototype,indices(rate_prototype))
  k5 = similar(rate_prototype,indices(rate_prototype)); k6 = similar(rate_prototype,indices(rate_prototype)); k7 = similar(rate_prototype,indices(rate_prototype)); k8 = similar(rate_prototype,indices(rate_prototype))
  k9 = similar(rate_prototype,indices(rate_prototype)); k10 = similar(rate_prototype,indices(rate_prototype)); k11 = similar(rate_prototype,indices(rate_prototype)); k12 = similar(rate_prototype,indices(rate_prototype))
  k13 = similar(rate_prototype,indices(rate_prototype)); update = similar(u); utilde = similar(u); k = similar(rate_prototype,indices(rate_prototype))
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits);
  TsitPap8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,update,tmp,atmp,k,tab)
end

alg_cache(alg::TsitPap8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = TsitPap8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable ExplicitRKCache{uType,rateType,uEltypeNoUnits,ksEltype,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  tmp::uType
  utilde::rateType
  uEEst::rateType
  atmp::uEltypeNoUnits
  uprev::uType
  fsalfirst::ksEltype
  fsallast::ksEltype
  kk::Vector{ksEltype}
  tab::TabType
end

function alg_cache(alg::ExplicitRK,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  kk = Vector{typeof(rate_prototype)}(0)
  for i = 1:alg.tableau.stages
    push!(kk,similar(rate_prototype,indices(rate_prototype)))
  end
  fsalfirst = kk[1]
  if isfsal(alg.tableau)
    fsallast = kk[end]
  else
    fsallast = similar(rate_prototype,indices(rate_prototype))
  end
  utilde = similar(rate_prototype,indices(rate_prototype))
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  uEEst = similar(rate_prototype,indices(rate_prototype))
  tab = ExplicitRKConstantCache(alg.tableau,rate_prototype)
  ExplicitRKCache(u,tmp,utilde,uEEst,atmp,uprev,fsalfirst,fsallast,kk,tab)
end

immutable ExplicitRKConstantCache{MType,VType,KType} <: OrdinaryDiffEqConstantCache
  A::MType
  c::VType
  α::VType
  αEEst::VType
  stages::Int
  kk::KType
end

function ExplicitRKConstantCache(tableau,rate_prototype)
  @unpack A,c,α,αEEst,stages = tableau
  A = A' # Transpose A to column major looping
  kk = Array{typeof(rate_prototype)}(stages) # Not ks since that's for integrator.opts.dense
  ExplicitRKConstantCache(A,c,α,αEEst,stages,kk)
end

alg_cache(alg::ExplicitRK,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = ExplicitRKConstantCache(alg.tableau,rate_prototype)

immutable Feagin10Cache{uType,uEltypeNoUnits,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  k17::rateType
  tmp::uType
  atmp::uEltypeNoUnits
  uprev::uType
  k::rateType
  tab::TabType
end

function alg_cache(alg::Feagin10,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = Feagin10ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype)); k2 = similar(rate_prototype,indices(rate_prototype)); k3 = similar(rate_prototype,indices(rate_prototype)); k4 = similar(rate_prototype,indices(rate_prototype)); k5 = similar(rate_prototype,indices(rate_prototype))
  k6 = similar(rate_prototype,indices(rate_prototype)); k7 = similar(rate_prototype,indices(rate_prototype)); k8 = similar(rate_prototype,indices(rate_prototype)); k9 = similar(rate_prototype,indices(rate_prototype)); k10 = similar(rate_prototype,indices(rate_prototype))
  k11 = similar(rate_prototype,indices(rate_prototype)); k12 = similar(rate_prototype,indices(rate_prototype)); k13 = similar(rate_prototype,indices(rate_prototype)); k14 = similar(rate_prototype,indices(rate_prototype))
  k15 = similar(rate_prototype,indices(rate_prototype)); k16 = similar(rate_prototype,indices(rate_prototype)); k17 = similar(rate_prototype,indices(rate_prototype))
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype,indices(rate_prototype))

  Feagin10Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,tmp,atmp,uprev,k,tab)
end

alg_cache(alg::Feagin10,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = Feagin10ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Feagin12Cache{uType,uEltypeNoUnits,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  k17::rateType
  k18::rateType
  k19::rateType
  k20::rateType
  k21::rateType
  k22::rateType
  k23::rateType
  k24::rateType
  k25::rateType
  tmp::uType
  atmp::uEltypeNoUnits
  uprev::uType
  k::rateType
  tab::TabType
end

function alg_cache(alg::Feagin12,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = Feagin12ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype)); k2 = similar(rate_prototype,indices(rate_prototype)); k3 = similar(rate_prototype,indices(rate_prototype)); k4 = similar(rate_prototype,indices(rate_prototype)); k5 = similar(rate_prototype,indices(rate_prototype))
  k6 = similar(rate_prototype,indices(rate_prototype)); k7 = similar(rate_prototype,indices(rate_prototype)); k8 = similar(rate_prototype,indices(rate_prototype)); k9 = similar(rate_prototype,indices(rate_prototype)); k10 = similar(rate_prototype,indices(rate_prototype))
  k11 = similar(rate_prototype,indices(rate_prototype)); k12 = similar(rate_prototype,indices(rate_prototype)); k13 = similar(rate_prototype,indices(rate_prototype)); k14 = similar(rate_prototype,indices(rate_prototype))
  k15 = similar(rate_prototype,indices(rate_prototype)); k16 = similar(rate_prototype,indices(rate_prototype)); k17 = similar(rate_prototype,indices(rate_prototype)); k18 = similar(rate_prototype,indices(rate_prototype))
  k19 = similar(rate_prototype,indices(rate_prototype)); k20 = similar(rate_prototype,indices(rate_prototype)); k21 = similar(rate_prototype,indices(rate_prototype)); k22 = similar(rate_prototype,indices(rate_prototype))
  k23 = similar(rate_prototype,indices(rate_prototype)); k24 = similar(rate_prototype,indices(rate_prototype)); k25 = similar(rate_prototype,indices(rate_prototype))
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype,indices(rate_prototype))

  Feagin12Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,tmp,atmp,uprev,k,tab)
end

alg_cache(alg::Feagin12,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = Feagin12ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))


immutable Feagin14Cache{uType,uEltypeNoUnits,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  k17::rateType
  k18::rateType
  k19::rateType
  k20::rateType
  k21::rateType
  k22::rateType
  k23::rateType
  k24::rateType
  k25::rateType
  k26::rateType
  k27::rateType
  k28::rateType
  k29::rateType
  k30::rateType
  k31::rateType
  k32::rateType
  k33::rateType
  k34::rateType
  k35::rateType
  tmp::uType
  atmp::uEltypeNoUnits
  uprev::uType
  k::rateType
  tab::TabType
end

function alg_cache(alg::Feagin14,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  tab = Feagin14ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype,indices(rate_prototype)); k2 = similar(rate_prototype,indices(rate_prototype)); k3 = similar(rate_prototype,indices(rate_prototype)); k4 = similar(rate_prototype,indices(rate_prototype)); k5 = similar(rate_prototype,indices(rate_prototype))
  k6 = similar(rate_prototype,indices(rate_prototype)); k7 = similar(rate_prototype,indices(rate_prototype)); k8 = similar(rate_prototype,indices(rate_prototype)); k9 = similar(rate_prototype,indices(rate_prototype)); k10 = similar(rate_prototype,indices(rate_prototype))
  k11 = similar(rate_prototype,indices(rate_prototype)); k12 = similar(rate_prototype,indices(rate_prototype)); k13 = similar(rate_prototype,indices(rate_prototype)); k14 = similar(rate_prototype,indices(rate_prototype))
  k15 = similar(rate_prototype,indices(rate_prototype)); k16 = similar(rate_prototype,indices(rate_prototype)); k17 = similar(rate_prototype,indices(rate_prototype)); k18 = similar(rate_prototype,indices(rate_prototype))
  k19 = similar(rate_prototype,indices(rate_prototype)); k20 = similar(rate_prototype,indices(rate_prototype)); k21 = similar(rate_prototype,indices(rate_prototype)); k22 = similar(rate_prototype,indices(rate_prototype))
  k23 = similar(rate_prototype,indices(rate_prototype)); k24 = similar(rate_prototype,indices(rate_prototype)); k25 = similar(rate_prototype,indices(rate_prototype))
  k26 = similar(rate_prototype,indices(rate_prototype)); k27 = similar(rate_prototype,indices(rate_prototype)); k28 = similar(rate_prototype,indices(rate_prototype))
  k29 = similar(rate_prototype,indices(rate_prototype)); k30 = similar(rate_prototype,indices(rate_prototype)); k31 = similar(rate_prototype,indices(rate_prototype)); k32 = similar(rate_prototype,indices(rate_prototype))
  k33 = similar(rate_prototype,indices(rate_prototype)); k34 = similar(rate_prototype,indices(rate_prototype)); k35 = similar(rate_prototype,indices(rate_prototype))
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype,indices(rate_prototype))

  Feagin14Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,
                k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,
                k31,k32,k33,k34,k35,tmp,atmp,uprev,k,tab)
end

alg_cache(alg::Feagin14,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = Feagin14ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))


immutable Rosenbrock23Cache{uType,rateType,vecuType,JType,TabType,TFType,UFType} <: OrdinaryDiffEqMutableCache
  u::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du1::rateType
  du2::rateType
  f₁::rateType
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uType
  J::JType
  W::JType
  tmp::uType
  tmp2::uType
  tab::TabType
  tf::TFType
  uf::UFType
end

immutable Rosenbrock32Cache{uType,rateType,vecuType,JType,TabType,TFType,UFType} <: OrdinaryDiffEqMutableCache
  u::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du1::rateType
  du2::rateType
  f₁::rateType
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uType
  J::JType
  W::JType
  tmp::uType
  tmp2::uType
  tab::TabType
  tf::TFType
  uf::UFType
end

function alg_cache(alg::Rosenbrock23,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  k₁ = similar(rate_prototype,indices(rate_prototype))
  k₂ = similar(rate_prototype,indices(rate_prototype))
  k₃ = similar(rate_prototype,indices(rate_prototype))
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  # f₀ = similar(u) fsalfirst
  f₁ = similar(rate_prototype,indices(rate_prototype))
  vectmp = similar(vec(u))
  vectmp2 = similar(vec(u))
  vectmp3 = similar(vec(u))
  fsalfirst = similar(rate_prototype,indices(rate_prototype))
  fsallast = similar(rate_prototype,indices(rate_prototype))
  dT = similar(u)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J); tmp2 = similar(u)
  tmp = reshape(vectmp2,size(u)...)
  tab = Rosenbrock23ConstantCache(uEltypeNoUnits,identity,identity)
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev,du2)
  uf = UJacobianWrapper(vfr,t)
  Rosenbrock23Cache(u,k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,tmp2,tab,tf,uf)
end

function alg_cache(alg::Rosenbrock32,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  k₁ = similar(rate_prototype,indices(rate_prototype))
  k₂ = similar(rate_prototype,indices(rate_prototype))
  k₃ = similar(rate_prototype,indices(rate_prototype))
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  # f₀ = similar(u) fsalfirst
  f₁ = similar(rate_prototype,indices(rate_prototype))
  vectmp = similar(vec(u))
  vectmp2 = similar(vec(u))
  vectmp3 = similar(vec(u))
  fsalfirst = similar(rate_prototype,indices(rate_prototype))
  fsallast = similar(rate_prototype,indices(rate_prototype))
  dT = similar(u)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J); tmp2 = similar(u)
  tmp = reshape(vectmp2,size(u)...)
  tab = Rosenbrock32ConstantCache(uEltypeNoUnits,identity,identity)
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev,du2)
  uf = UJacobianWrapper(vfr,t)
  Rosenbrock32Cache(u,k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,tmp2,tab,tf,uf)
end

immutable Rosenbrock23ConstantCache{T,TF,UF} <: OrdinaryDiffEqConstantCache
  c₃₂::T
  d::T
  tf::TF
  uf::UF
end

function Rosenbrock23ConstantCache(T::Type,tf,uf)
  c₃₂ = T(6 + sqrt(2))
  d = T(1/(2+sqrt(2)))
  Rosenbrock23ConstantCache(c₃₂,d,tf,uf)
end

function alg_cache(alg::Rosenbrock23,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock23ConstantCache(uEltypeNoUnits,tf,uf)
end

immutable Rosenbrock32ConstantCache{T,TF,UF} <: OrdinaryDiffEqConstantCache
  c₃₂::T
  d::T
  tf::TF
  uf::UF
end

function Rosenbrock32ConstantCache(T::Type,tf,uf)
  c₃₂ = T(6 + sqrt(2))
  d = T(1/(2+sqrt(2)))
  Rosenbrock32ConstantCache(c₃₂,d,tf,uf)
end

function alg_cache(alg::Rosenbrock32,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock32ConstantCache(uEltypeNoUnits,tf,uf)
end

immutable ImplicitEulerCache{uType,vecuType,DiffCacheType,rateType,rhsType,adfType,CS} <: OrdinaryDiffEqMutableCache
  u::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  u_old::uType
  k::rateType
  fsalfirst::rateType
  rhs::rhsType
  adf::adfType
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  u_old = similar(u); k = similar(rate_prototype,indices(rate_prototype))
  dual_cache = DiffCache(u,Val{determine_chunksize(u,alg)})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  rhs = RHS_IE(f,u_old,t,t,dual_cache,size(u),eachindex(u))
  fsalfirst = similar(rate_prototype,indices(rate_prototype))
  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  else
    adf = nothing
  end
  ImplicitEulerCache{typeof(u),typeof(uhold),typeof(dual_cache),typeof(k),typeof(rhs),typeof(adf),determine_chunksize(u,alg)}(u,uhold,dual_cache,u_old,k,fsalfirst,rhs,adf)
end

immutable ImplicitEulerConstantCache{vecuType,rhsType,adfType,CS} <: OrdinaryDiffEqMutableCache
  uhold::vecuType
  u_old::vecuType
  rhs::rhsType
  adf::adfType
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  u_old = Vector{typeof(u)}(1)
  rhs = RHS_IE_Scalar(f,u_old,t,t)
  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  else
    adf = nothing
  end
  ImplicitEulerConstantCache{typeof(uhold),typeof(rhs),typeof(adf),1}(uhold,u_old,rhs,adf)
end

immutable TrapezoidCache{uType,vecuType,DiffCacheType,rateType,rhsType,adfType,CS} <: OrdinaryDiffEqMutableCache
  u::uType
  uhold::vecuType
  u_old::uType
  f_old::rateType
  dual_cache::DiffCacheType
  k::rateType
  rhs::rhsType
  adf::adfType
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  u_old = similar(u); k = similar(rate_prototype,indices(rate_prototype))
  uhold = vec(u); f_old = similar(rate_prototype,indices(rate_prototype))
  dual_cache = DiffCache(u,Val{determine_chunksize(u,alg)})
  rhs = RHS_Trap(f,u_old,f_old,t,t,size(u),dual_cache,eachindex(u))
  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  else
    adf = nothing
  end
  TrapezoidCache{typeof(u),typeof(uhold),typeof(dual_cache),typeof(k),
    typeof(rhs),typeof(adf),determine_chunksize(u,alg)}(
    u,uhold,u_old,f_old,dual_cache,k,rhs,adf)
end


immutable TrapezoidConstantCache{vecuType,rhsType,adfType,CS} <: OrdinaryDiffEqMutableCache
  uhold::vecuType
  u_old::vecuType
  rhs::rhsType
  adf::adfType
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  u_old = Vector{typeof(u)}(1)
  rhs = RHS_Trap_Scalar(f,u_old,rate_prototype,t,t)
  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  else
    adf = nothing
  end
  TrapezoidConstantCache{typeof(uhold),typeof(rhs),typeof(adf),1}(uhold,u_old,rhs,adf)
end



get_chunksize(cache::DECache) = error("This cache does not have a chunksize.")
get_chunksize{uType,DiffCacheType,rateType,CS}(cache::ImplicitEulerCache{uType,DiffCacheType,rateType,CS}) = CS
get_chunksize{uType,DiffCacheType,rateType,CS}(cache::TrapezoidCache{uType,DiffCacheType,rateType,CS}) = CS
get_chunksize{CS}(cache::ODEChunkCache{CS}) = CS
