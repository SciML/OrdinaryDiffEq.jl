@cache struct TanYam7Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
  atmp::uNoUnitsType
  k::rateType
  tab::TabType
end

function alg_cache(alg::TanYam7,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = TanYam7ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zero(rate_prototype); k2 = zero(rate_prototype) ; k3 = zero(rate_prototype); k4 = zero(rate_prototype)
  k5 = zero(rate_prototype); k6 = zero(rate_prototype) ; k7 = zero(rate_prototype); k8 = zero(rate_prototype)
  k9 = zero(rate_prototype); k10= zero(rate_prototype) ;
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = zero(rate_prototype)
  TanYam7Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,tmp,atmp,k,tab)
end

alg_cache(alg::TanYam7,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = TanYam7ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

@cache struct DP8Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
  kupdate::rateType
  udiff::rateType
  bspl::rateType
  dense_tmp3::rateType
  dense_tmp4::rateType
  dense_tmp5::rateType
  dense_tmp6::rateType
  dense_tmp7::rateType
  utilde::uType
  tmp::uType
  atmp::uNoUnitsType
  tab::TabType
end

function alg_cache(alg::DP8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  k1 = zero(rate_prototype); k2  = zero(rate_prototype); k3  = zero(rate_prototype);  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype); k6  = zero(rate_prototype); k7  = zero(rate_prototype);  k8 = zero(rate_prototype)
  k9 = zero(rate_prototype); k10 = zero(rate_prototype); k11 = zero(rate_prototype); k12 = zero(rate_prototype)
  kupdate = zero(rate_prototype); utilde = similar(u);
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  k13 = zero(rate_prototype)
  k14 = zero(rate_prototype)
  k15 = zero(rate_prototype)
  k16 = zero(rate_prototype)
  udiff = zero(rate_prototype)
  bspl = zero(rate_prototype)
  # dense_tmp1 = udiff
  # dense_tmp2 = bspl
  dense_tmp3 = zero(rate_prototype)
  dense_tmp4 = zero(rate_prototype)
  dense_tmp5 = zero(rate_prototype)
  dense_tmp6 = zero(rate_prototype)
  dense_tmp7 = zero(rate_prototype)
  tab = DP8ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  DP8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,kupdate,
           udiff,bspl,dense_tmp3,dense_tmp4,dense_tmp5,dense_tmp6,dense_tmp7,
           utilde,tmp,atmp,tab)
end

alg_cache(alg::DP8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = DP8ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

@cache struct TsitPap8Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
  tmp::uType
  atmp::uNoUnitsType
  k::rateType
  tab::TabType
end

function alg_cache(alg::TsitPap8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = TsitPap8ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype); k4 = zero(rate_prototype)
  k5 = zero(rate_prototype); k6 = zero(rate_prototype); k7 = zero(rate_prototype); k8 = zero(rate_prototype)
  k9 = zero(rate_prototype); k10 = zero(rate_prototype); k11 = zero(rate_prototype); k12 = zero(rate_prototype)
  k13 = zero(rate_prototype); utilde = similar(u); k = zero(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits);
  TsitPap8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,tmp,atmp,k,tab)
end

alg_cache(alg::TsitPap8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = TsitPap8ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
