
################################################################################

# Shampine's Low-order Rosenbrocks

type Rosenbrock23Cache{uType,uArrayType,rateType,du2Type,LinuType,vecuType,JType,TabType,TFType,UFType,F,JCType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du1::rateType
  du2::du2Type
  f₁::rateType
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uArrayType
  J::JType
  W::JType
  tmp::uArrayType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::LinuType
  linsolve_tmp_vec::vecuType
  linsolve::F
  jac_config::JCType
end

u_cache(c::Rosenbrock23Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock23Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock23Cache) = (c.J,c.W)
vecu_cache(c::Rosenbrock23Cache) = (c.vectmp,c.vectmp2,c.vectmp3)

type Rosenbrock32Cache{uType,uArrayType,rateType,du2Type,LinuType,vecuType,JType,TabType,TFType,UFType,F,JCType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du1::rateType
  du2::du2Type
  f₁::rateType
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uArrayType
  J::JType
  W::JType
  tmp::uArrayType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::LinuType
  linsolve_tmp_vec::vecuType
  linsolve::F
  jac_config::JCType
end

u_cache(c::Rosenbrock32Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock32Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock32Cache) = (c.J,c.W)
vecu_cache(c::Rosenbrock32Cache) = (c.vectmp,c.vectmp2,c.vectmp3)

function alg_cache(alg::Rosenbrock23,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  # f₀ = similar(u) fsalfirst
  f₁ = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J);
  tmp = similar(u,indices(u))
  tab = Rosenbrock23ConstantCache(uEltypeNoUnits,identity,identity)
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  Rosenbrock23Cache(u,uprev,k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,
                    fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,linsolve_tmp_vec,
                    alg.linsolve,jac_config)
end

function alg_cache(alg::Rosenbrock32,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  # f₀ = similar(u) fsalfirst
  f₁ = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J); tmp = similar(u,indices(u))
  tab = Rosenbrock32ConstantCache(uEltypeNoUnits,identity,identity)
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  Rosenbrock32Cache(u,uprev,k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,linsolve_tmp_vec,alg.linsolve,jac_config)
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

function alg_cache(alg::Rosenbrock23,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
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

function alg_cache(alg::Rosenbrock32,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock32ConstantCache(uEltypeNoUnits,tf,uf)
end

################################################################################

### 3rd order specialized Rosenbrocks

immutable Rosenbrock33ConstantCache{TF,UF,Tab} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
end

type Rosenbrock33Cache{uType,uArrayType,rateType,du2Type,LinuType,vecuType,JType,TabType,TFType,UFType,F,JCType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du::rateType
  du1::rateType
  du2::du2Type
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  vectmp4::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uArrayType
  J::JType
  W::JType
  tmp::uArrayType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::LinuType
  linsolve_tmp_vec::vecuType
  linsolve::F
  jac_config::JCType
end

u_cache(c::Rosenbrock33Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock33Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock33Cache) = (c.J,c.W)
vecu_cache(c::Rosenbrock33Cache) = (c.vectmp,c.vectmp2,c.vectmp3)

function alg_cache(alg::ROS3P,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  du = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  vectmp4 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J);
  tmp = similar(u,indices(u))
  tab = ROS3PConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  Rosenbrock33Cache(u,uprev,du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve_tmp_vec,alg.linsolve,jac_config)
end

function alg_cache(alg::ROS3P,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock33ConstantCache(tf,uf,ROS3PConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits)))
end

type Rosenbrock34Cache{uType,uArrayType,rateType,du2Type,LinuType,vecuType,JType,TabType,TFType,UFType,F,JCType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du::rateType
  du1::rateType
  du2::du2Type
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  vectmp4::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uArrayType
  J::JType
  W::JType
  tmp::uArrayType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::LinuType
  linsolve_tmp_vec::vecuType
  linsolve::F
  jac_config::JCType
end

u_cache(c::Rosenbrock34Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock34Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock34Cache) = (c.J,c.W)
vecu_cache(c::Rosenbrock34Cache) = (c.vectmp,c.vectmp2,c.vectmp3)

function alg_cache(alg::Rodas3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  du = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  vectmp4 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J);
  tmp = similar(u,indices(u))
  tab = Rodas3ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  Rosenbrock34Cache(u,uprev,du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve_tmp_vec,alg.linsolve,jac_config)
end

immutable Rosenbrock34ConstantCache{TF,UF,Tab} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
end

function alg_cache(alg::Rodas3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock34ConstantCache(tf,uf,Rodas3ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits)))
end


################################################################################

immutable Rosenbrock4ConstantCache{TF,UF,Tab} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
end

type Rosenbrock4Cache{uType,uArrayType,rateType,du2Type,LinuType,vecuType,JType,TabType,TFType,UFType,F,JCType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du::rateType
  du1::rateType
  du2::du2Type
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  vectmp4::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uArrayType
  J::JType
  W::JType
  tmp::uArrayType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::LinuType
  linsolve_tmp_vec::vecuType
  linsolve::F
  jac_config::JCType
end

u_cache(c::Rosenbrock4Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock4Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock4Cache) = (c.J,c.W)
vecu_cache(c::Rosenbrock4Cache) = (c.vectmp,c.vectmp2,c.vectmp3)

function alg_cache(alg::RosShamp4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  du = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  vectmp4 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J);
  tmp = similar(u,indices(u))
  tab = RosShamp4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  Rosenbrock4Cache(u,uprev,du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve_tmp_vec,alg.linsolve,jac_config)
end

function alg_cache(alg::RosShamp4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock4ConstantCache(tf,uf,RosShamp4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits)))
end

function alg_cache(alg::Veldd4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  du = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  vectmp4 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J);
  tmp = similar(u,indices(u))
  tab = Veldd4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  Rosenbrock4Cache(u,uprev,du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve_tmp_vec,alg.linsolve,jac_config)
end

function alg_cache(alg::Veldd4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock4ConstantCache(tf,uf,Veldd4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits)))
end

function alg_cache(alg::Velds4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  du = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  vectmp4 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J);
  tmp = similar(u,indices(u))
  tab = Velds4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  Rosenbrock4Cache(u,uprev,du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve_tmp_vec,alg.linsolve,jac_config)
end

function alg_cache(alg::Velds4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock4ConstantCache(tf,uf,Velds4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits)))
end

function alg_cache(alg::GRK4T,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  du = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  vectmp4 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J);
  tmp = similar(u,indices(u))
  tab = GRK4TConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  Rosenbrock4Cache(u,uprev,du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve_tmp_vec,alg.linsolve,jac_config)
end

function alg_cache(alg::GRK4T,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock4ConstantCache(tf,uf,GRK4TConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits)))
end

function alg_cache(alg::GRK4A,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  du = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  vectmp4 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J);
  tmp = similar(u,indices(u))
  tab = GRK4AConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  Rosenbrock4Cache(u,uprev,du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve_tmp_vec,alg.linsolve,jac_config)
end

function alg_cache(alg::GRK4A,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock4ConstantCache(tf,uf,GRK4AConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits)))
end

function alg_cache(alg::Ros4LStab,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  du = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  vectmp4 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J);
  tmp = similar(u,indices(u))
  tab = Ros4LStabConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  Rosenbrock4Cache(u,uprev,du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve_tmp_vec,alg.linsolve,jac_config)
end

function alg_cache(alg::Ros4LStab,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock4ConstantCache(tf,uf,Ros4LStabConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits)))
end

###############################################################################

### Rodas methods

immutable Rodas4ConstantCache{TF,UF,Tab} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
end

type Rodas4Cache{uType,uArrayType,rateType,du2Type,LinuType,vecuType,JType,TabType,TFType,UFType,F,JCType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du::rateType
  du1::rateType
  du2::du2Type
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  vectmp4::vecuType
  vectmp5::vecuType
  vectmp6::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uArrayType
  J::JType
  W::JType
  tmp::uArrayType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::LinuType
  linsolve_tmp_vec::vecuType
  linsolve::F
  jac_config::JCType
end

u_cache(c::Rodas4Cache) = (c.dT,c.tmp)
du_cache(c::Rodas4Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rodas4Cache) = (c.J,c.W)
vecu_cache(c::Rodas4Cache) = (c.vectmp,c.vectmp2,c.vectmp3)

function alg_cache(alg::Rodas4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  du = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  vectmp4 = vec(similar(u,indices(u)))
  vectmp5 = vec(similar(u,indices(u)))
  vectmp6 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J);
  tmp = similar(u,indices(u))
  tab = Rodas4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  Rodas4Cache(u,uprev,k1,k2,k3,k4,du,du1,du2,vectmp,vectmp2,vectmp3,vectmp4,
                    vectmp5,vectmp6,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve_tmp_vec,alg.linsolve,jac_config)
end

function alg_cache(alg::Rodas4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rodas4ConstantCache(tf,uf,Rodas4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits)))
end

function alg_cache(alg::Rodas42,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rodas4ConstantCache(tf,uf,Rodas42ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits)))
end

function alg_cache(alg::Rodas4P,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rodas4ConstantCache(tf,uf,Rodas4PConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits)))
end
