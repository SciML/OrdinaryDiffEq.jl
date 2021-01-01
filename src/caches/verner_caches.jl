@cache struct Vern6Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uType
  tmp::uType
  rtmp::rateType
  atmp::uNoUnitsType
  tab::TabType
end

function alg_cache(alg::Vern6,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = Vern6Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype); k3 = k2;                   k4 = zero(rate_prototype);
  k5 = zero(rate_prototype); k6 = zero(rate_prototype); k7 = zero(rate_prototype);
  k8 = k3;                   k9 = zero(rate_prototype);
  utilde = zero(u); tmp = zero(u); atmp = similar(u,uEltypeNoUnits);
  rtmp = zero(rate_prototype);
  Vern6Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,utilde,tmp,rtmp,atmp,tab)
end

struct Vern6ConstantCache{TabType} <: OrdinaryDiffEqConstantCache
  tab::TabType
end

function alg_cache(alg::Vern6,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = Vern6Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  Vern6ConstantCache(tab)
end

@cache struct Vern7Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uType
  tmp::uType
  rtmp::rateType
  atmp::uNoUnitsType
  tab::TabType
end

function alg_cache(alg::Vern7,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = Vern7Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = k2                  ; k4 = zero(rate_prototype);
  k5 = zero(rate_prototype); k6 = zero(rate_prototype); k7 = zero(rate_prototype); k8 = zero(rate_prototype);
  k9 = zero(rate_prototype); k10 = k2                 ; utilde = zero(u)
  tmp = zero(u); atmp = similar(u,uEltypeNoUnits)
  rtmp = zero(rate_prototype)
  Vern7Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,tmp,rtmp,atmp,tab)
end

struct Vern7ConstantCache{TabType} <: OrdinaryDiffEqConstantCache
  tab::TabType
end

function alg_cache(alg::Vern7,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = Vern7Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  Vern7ConstantCache(tab)
end

@cache struct Vern8Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uType
  tmp::uType
  rtmp::rateType
  atmp::uNoUnitsType
  tab::TabType
end

function alg_cache(alg::Vern8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = Vern8Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = k2;
  k4 = zero(rate_prototype);
  k5 = k2                  ; k6 = zero(rate_prototype); k7 = zero(rate_prototype);
  k8 = zero(rate_prototype); tmp = zero(u)
  k9 = zero(rate_prototype); k10 = zero(rate_prototype); k11 = zero(rate_prototype);
  k12 = zero(rate_prototype); k13 = k4
  utilde = zero(u)
  atmp = similar(u,uEltypeNoUnits)
  rtmp = zero(rate_prototype)
  Vern8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,tmp,rtmp,atmp,tab)
end

struct Vern8ConstantCache{TabType} <: OrdinaryDiffEqConstantCache
  tab::TabType
end

function alg_cache(alg::Vern8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = Vern8Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  Vern8ConstantCache(tab)
end

@cache struct Vern9Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uType
  tmp::uType
  rtmp::rateType
  atmp::uNoUnitsType
  tab::TabType
end

function alg_cache(alg::Vern9,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = Vern9Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  k1 = zero(rate_prototype); k2 = zero(rate_prototype);k3 = k2;
  k4 = zero(rate_prototype);
  k5 = k3;                   k6 = zero(rate_prototype);k7 = k4;
  k8 = k5;
  k9 = zero(rate_prototype); k10 = zero(rate_prototype); k11 = zero(rate_prototype);
  k12 = zero(rate_prototype);
  k13 = zero(rate_prototype); k14 = zero(rate_prototype); k15 = zero(rate_prototype);
  k16 =k6;
  utilde = zero(u); tmp = zero(u); atmp = similar(u,uEltypeNoUnits);
  rtmp = zero(rate_prototype)
  Vern9Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,utilde,tmp,rtmp,atmp,tab)
end

struct Vern9ConstantCache{TabType} <: OrdinaryDiffEqConstantCache
  tab::TabType
end

function alg_cache(alg::Vern9,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = Vern9Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  Vern9ConstantCache(tab)
end