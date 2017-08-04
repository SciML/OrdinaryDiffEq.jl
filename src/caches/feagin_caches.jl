struct Feagin10Cache{uType,uEltypeNoUnits,rateType,TabType} <: OrdinaryDiffEqMutableCache
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
  k14::rateType
  k15::rateType
  k16::rateType
  k17::rateType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

u_cache(c::Feagin10Cache) = (c.atmp,)
du_cache(c::Feagin10Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.k17,c.k)

function alg_cache(alg::Feagin10,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})
  tab = Feagin10ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype); k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype); k7 = zeros(rate_prototype); k8 = zeros(rate_prototype); k9 = zeros(rate_prototype); k10 = zeros(rate_prototype)
  k11 = zeros(rate_prototype); k12 = zeros(rate_prototype); k13 = zeros(rate_prototype); k14 = zeros(rate_prototype)
  k15 = zeros(rate_prototype); k16 = zeros(rate_prototype); k17 = zeros(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u)); k = zeros(rate_prototype)

  Feagin10Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,tmp,atmp,k,tab)
end

alg_cache(alg::Feagin10,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{false}}) = Feagin10ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))

struct Feagin12Cache{uType,uEltypeNoUnits,rateType,TabType} <: OrdinaryDiffEqMutableCache
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
  k::rateType
  tab::TabType
end

u_cache(c::Feagin12Cache) = (c.atmp,)
du_cache(c::Feagin12Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.k17,c.k18,c.k19,c.k20,c.k21,c.k22,c.k23,c.k24,c.k25,c.k)

function alg_cache(alg::Feagin12,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})
  tab = Feagin12ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype); k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype); k7 = zeros(rate_prototype); k8 = zeros(rate_prototype); k9 = zeros(rate_prototype); k10 = zeros(rate_prototype)
  k11 = zeros(rate_prototype); k12 = zeros(rate_prototype); k13 = zeros(rate_prototype); k14 = zeros(rate_prototype)
  k15 = zeros(rate_prototype); k16 = zeros(rate_prototype); k17 = zeros(rate_prototype); k18 = zeros(rate_prototype)
  k19 = zeros(rate_prototype); k20 = zeros(rate_prototype); k21 = zeros(rate_prototype); k22 = zeros(rate_prototype)
  k23 = zeros(rate_prototype); k24 = zeros(rate_prototype); k25 = zeros(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u)); k = zeros(rate_prototype)

  Feagin12Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,tmp,atmp,k,tab)
end

alg_cache(alg::Feagin12,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{false}}) = Feagin12ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))


struct Feagin14Cache{uType,uEltypeNoUnits,rateType,TabType} <: OrdinaryDiffEqMutableCache
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
  k::rateType
  tab::TabType
end

u_cache(c::Feagin14Cache) = (c.atmp,)
du_cache(c::Feagin14Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.k17,c.k18,c.k19,c.k20,c.k21,c.k22,c.k23,c.k24,c.k25,c.k26,c.k27,c.k28,c.k29,c.k30,c.k31,c.k32,c.k33,c.k34,c.k35,c.k)


function alg_cache(alg::Feagin14,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})
  tab = Feagin14ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype); k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype); k7 = zeros(rate_prototype); k8 = zeros(rate_prototype); k9 = zeros(rate_prototype); k10 = zeros(rate_prototype)
  k11 = zeros(rate_prototype); k12 = zeros(rate_prototype); k13 = zeros(rate_prototype); k14 = zeros(rate_prototype)
  k15 = zeros(rate_prototype); k16 = zeros(rate_prototype); k17 = zeros(rate_prototype); k18 = zeros(rate_prototype)
  k19 = zeros(rate_prototype); k20 = zeros(rate_prototype); k21 = zeros(rate_prototype); k22 = zeros(rate_prototype)
  k23 = zeros(rate_prototype); k24 = zeros(rate_prototype); k25 = zeros(rate_prototype)
  k26 = zeros(rate_prototype); k27 = zeros(rate_prototype); k28 = zeros(rate_prototype)
  k29 = zeros(rate_prototype); k30 = zeros(rate_prototype); k31 = zeros(rate_prototype); k32 = zeros(rate_prototype)
  k33 = zeros(rate_prototype); k34 = zeros(rate_prototype); k35 = zeros(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u)); k = zeros(rate_prototype)

  Feagin14Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,
                k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,
                k31,k32,k33,k34,k35,tmp,atmp,k,tab)
end

alg_cache(alg::Feagin14,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{false}}) = Feagin14ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))
