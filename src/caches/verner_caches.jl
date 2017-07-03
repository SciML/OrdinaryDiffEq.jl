immutable Vern6Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern6Cache) = (c.utilde,c.atmp)
du_cache(c::Vern6Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9)

function alg_cache(alg::Vern6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern6ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype);
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype); k7 = zeros(rate_prototype);
  k8 = zeros(rate_prototype); k9 = zeros(rate_prototype);
  utilde = similar(u,indices(u)); tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u));
  Vern6Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,utilde,tmp,atmp,tab)
end

alg_cache(alg::Vern6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern6ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Vern7Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
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
  utilde::uArrayType
  update::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern7Cache) = (c.utilde,c.update,c.atmp)
du_cache(c::Vern7Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10)

function alg_cache(alg::Vern7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype);
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype); k7 = zeros(rate_prototype); k8 = zeros(rate_prototype);
  k9 = zeros(rate_prototype); k10 = zeros(rate_prototype); utilde = similar(u,indices(u)); update = similar(u,indices(u))
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u))
  Vern7Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))


immutable Vern8Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
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
  utilde::uArrayType
  update::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern8Cache) = (c.utilde,c.update,c.atmp)
du_cache(c::Vern8Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13)

function alg_cache(alg::Vern8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype);
  k4 = zeros(rate_prototype);
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype); k7 = zeros(rate_prototype);
  k8 = zeros(rate_prototype); tmp = similar(u)
  k9 = zeros(rate_prototype); k10 = zeros(rate_prototype); k11 = zeros(rate_prototype);
  k12 = zeros(rate_prototype); k13 = zeros(rate_prototype)
  utilde = similar(u,indices(u)); update = similar(u,indices(u));
  atmp = similar(u,uEltypeNoUnits,indices(u))
  Vern8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Vern9Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
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
  utilde::uArrayType
  update::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern9Cache) = (c.utilde,c.update,c.atmp)
du_cache(c::Vern9Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16)

function alg_cache(alg::Vern9,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern9ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype);k3 = zeros(rate_prototype);
  k4 = zeros(rate_prototype);
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype);k7 = zeros(rate_prototype);
  k8 = zeros(rate_prototype);
  k9 = zeros(rate_prototype); k10 = zeros(rate_prototype); k11 = zeros(rate_prototype);
  k12 = zeros(rate_prototype); update = similar(u,indices(u))
  k13 = zeros(rate_prototype); k14 = zeros(rate_prototype); k15 = zeros(rate_prototype);
  k16 =zeros(rate_prototype);
  utilde = similar(u,indices(u)); tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u));
  Vern9Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern9,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern9ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
