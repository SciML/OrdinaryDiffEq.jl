
struct TanYam7Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

u_cache(c::TanYam7Cache) = (c.utilde,c.atmp)
du_cache(c::TanYam7Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k)

function alg_cache(alg::TanYam7,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tab = TanYam7ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype) ; k3 = zeros(rate_prototype); k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype) ; k7 = zeros(rate_prototype); k8 = zeros(rate_prototype)
  k9 = zeros(rate_prototype); k10= zeros(rate_prototype) ;
  utilde = similar(u,indices(u)); tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u)); k = zeros(rate_prototype)
  TanYam7Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,tmp,atmp,k,tab)
end

alg_cache(alg::TanYam7,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = TanYam7ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

struct DP8Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DP8Cache) = (c.utilde,c.atmp)
du_cache(c::DP8Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.kupdate,c.udiff,c.bspl,c.dense_tmp3,c.dense_tmp4,c.dense_tmp5,c.dense_tmp6,c.dense_tmp7)

function alg_cache(alg::DP8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  k1 = zeros(rate_prototype); k2  = zeros(rate_prototype); k3  = zeros(rate_prototype);  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype); k6  = zeros(rate_prototype); k7  = zeros(rate_prototype);  k8 = zeros(rate_prototype)
  k9 = zeros(rate_prototype); k10 = zeros(rate_prototype); k11 = zeros(rate_prototype); k12 = zeros(rate_prototype)
  kupdate = zeros(rate_prototype); utilde = similar(u,indices(u));
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u))
  k13 = zeros(rate_prototype)
  k14 = zeros(rate_prototype)
  k15 = zeros(rate_prototype)
  k16 = zeros(rate_prototype)
  udiff = zeros(rate_prototype)
  bspl = zeros(rate_prototype)
  # dense_tmp1 = udiff
  # dense_tmp2 = bspl
  dense_tmp3 = zeros(rate_prototype)
  dense_tmp4 = zeros(rate_prototype)
  dense_tmp5 = zeros(rate_prototype)
  dense_tmp6 = zeros(rate_prototype)
  dense_tmp7 = zeros(rate_prototype)
  tab = DP8ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  DP8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,kupdate,
           udiff,bspl,dense_tmp3,dense_tmp4,dense_tmp5,dense_tmp6,dense_tmp7,
           utilde,tmp,atmp,tab)
end

alg_cache(alg::DP8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = DP8ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

struct TsitPap8Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

u_cache(c::TsitPap8Cache) = (c.utilde,c.atmp)
du_cache(c::TsitPap8Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k)

function alg_cache(alg::TsitPap8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tab = TsitPap8ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype); k7 = zeros(rate_prototype); k8 = zeros(rate_prototype)
  k9 = zeros(rate_prototype); k10 = zeros(rate_prototype); k11 = zeros(rate_prototype); k12 = zeros(rate_prototype)
  k13 = zeros(rate_prototype); utilde = similar(u,indices(u)); k = zeros(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u));
  TsitPap8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,tmp,atmp,k,tab)
end

alg_cache(alg::TsitPap8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = TsitPap8ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
