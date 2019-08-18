abstract type RosenbrockMutableCache <: OrdinaryDiffEqMutableCache end
################################################################################

# Shampine's Low-order Rosenbrocks

@cache mutable struct Rosenbrock23Cache{uType,rateType,uNoUnitsType,N,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
  u::uType
  uprev::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du1::rateType
  du2::rateType
  f₁::rateType
  fsalfirst::rateType
  fsallast::rateType
  dT::rateType
  nlsolver::N
  tmp::rateType
  atmp::uNoUnitsType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

@cache mutable struct Rosenbrock32Cache{uType,rateType,uNoUnitsType,N,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
  u::uType
  uprev::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du1::rateType
  du2::rateType
  f₁::rateType
  fsalfirst::rateType
  fsallast::rateType
  dT::rateType
  nlsolver::N
  tmp::rateType
  atmp::uNoUnitsType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

function alg_cache(alg::Rosenbrock23,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  k₁ = zero(rate_prototype)
  k₂ = zero(rate_prototype)
  k₃ = zero(rate_prototype)
  du1 = zero(rate_prototype)
  du2 = zero(rate_prototype)
  # f₀ = similar(u) fsalfirst
  f₁ = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  fsallast = zero(rate_prototype)
  dT = zero(rate_prototype)
  J,W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  tmp = zero(rate_prototype)
  atmp = similar(u, uEltypeNoUnits)
  tab = Rosenbrock23Tableau(real(uBottomEltypeNoUnits))
  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)

  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  nlsolver = SemiImplicitNLSolver(W,J,du1,uf,jac_config)

  Rosenbrock23Cache(u,uprev,k₁,k₂,k₃,du1,du2,f₁,
                    fsalfirst,fsallast,dT,nlsolver,tmp,atmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Rosenbrock32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  k₁ = zero(rate_prototype)
  k₂ = zero(rate_prototype)
  k₃ = zero(rate_prototype)
  du1 = zero(rate_prototype)
  du2 = zero(rate_prototype)
  # f₀ = similar(u) fsalfirst
  f₁ = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  fsallast = zero(rate_prototype)
  dT = zero(rate_prototype)
  J,W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  tmp = zero(rate_prototype)
  atmp = similar(u, uEltypeNoUnits)
  tab = Rosenbrock32Tableau(real(uBottomEltypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  nlsolver = SemiImplicitNLSolver(W,J,du1,uf,jac_config)
  Rosenbrock32Cache(u,uprev,k₁,k₂,k₃,du1,du2,f₁,fsalfirst,fsallast,dT,nlsolver,tmp,atmp,tab,tf,uf,linsolve_tmp,linsolve,jac_config,grad_config)
end

struct Rosenbrock23ConstantCache{T,TF,UF,N,F} <: OrdinaryDiffEqConstantCache
  c₃₂::T
  d::T
  tf::TF
  uf::UF
  nlsolver::N
  linsolve::F
end

function Rosenbrock23ConstantCache(T::Type,tf,uf,nlsolver,linsolve)
  tab = Rosenbrock23Tableau(T)
  Rosenbrock23ConstantCache(tab.c₃₂,tab.d,tf,uf,nlsolver,linsolve)
end

function alg_cache(alg::Rosenbrock23,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  J,W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  linsolve = alg.linsolve(Val{:init},uf,u)
  nlsolver = SemiImplicitNLSolver(W,J,nothing,uf,nothing)
  Rosenbrock23ConstantCache(real(uBottomEltypeNoUnits),tf,uf,nlsolver,linsolve)
end

struct Rosenbrock32ConstantCache{T,TF,UF,N,F} <: OrdinaryDiffEqConstantCache
  c₃₂::T
  d::T
  tf::TF
  uf::UF
  nlsolver::N
  linsolve::F
end

function Rosenbrock32ConstantCache(T::Type,tf,uf,nlsolver,linsolve)
  tab = Rosenbrock32Tableau(T)
  Rosenbrock32ConstantCache(tab.c₃₂,tab.d,tf,uf,nlsolver,linsolve)
end

function alg_cache(alg::Rosenbrock32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  J,W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  linsolve = alg.linsolve(Val{:init},uf,u)
  nlsolver = SemiImplicitNLSolver(W,J,nothing,uf,nothing)
  Rosenbrock32ConstantCache(real(uBottomEltypeNoUnits),tf,uf,nlsolver,linsolve)
end

################################################################################

### 3rd order specialized Rosenbrocks

struct Rosenbrock33ConstantCache{TF,UF,Tab,N,F} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
  nlsolver::N
  linsolve::F
end

@cache mutable struct Rosenbrock33Cache{uType,rateType,uNoUnitsType,N,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
  u::uType
  uprev::uType
  du::rateType
  du1::rateType
  du2::rateType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  fsalfirst::rateType
  fsallast::rateType
  dT::rateType
  nlsolver::N
  tmp::rateType
  atmp::uNoUnitsType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

function alg_cache(alg::ROS3P,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  du = zero(rate_prototype)
  du1 = zero(rate_prototype)
  du2 = zero(rate_prototype)
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  fsallast = zero(rate_prototype)
  dT = zero(rate_prototype)
  J,W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  tmp = zero(rate_prototype)
  atmp = similar(u, uEltypeNoUnits)
  tab = ROS3PTableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  nlsolver = SemiImplicitNLSolver(W,J,du1,uf,jac_config)
  Rosenbrock33Cache(u,uprev,du,du1,du2,k1,k2,k3,k4,
                    fsalfirst,fsallast,dT,nlsolver,tmp,atmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::ROS3P,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  J,W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  linsolve = alg.linsolve(Val{:init},uf,u)
  nlsolver = SemiImplicitNLSolver(W,J,nothing,uf,nothing)
  Rosenbrock33ConstantCache(tf,uf,ROS3PTableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits)),nlsolver,linsolve)
end

@cache mutable struct Rosenbrock34Cache{uType,rateType,uNoUnitsType,N,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
  u::uType
  uprev::uType
  du::rateType
  du1::rateType
  du2::rateType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  fsalfirst::rateType
  fsallast::rateType
  dT::rateType
  nlsolver::N
  tmp::rateType
  atmp::uNoUnitsType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

function alg_cache(alg::Rodas3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  du = zero(rate_prototype)
  du1 = zero(rate_prototype)
  du2 = zero(rate_prototype)
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  fsallast = zero(rate_prototype)
  dT = zero(rate_prototype)
  J,W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  tmp = zero(rate_prototype)
  atmp = similar(u, uEltypeNoUnits)
  tab = Rodas3Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  nlsolver = SemiImplicitNLSolver(W,J,du1,uf,jac_config)
  Rosenbrock34Cache(u,uprev,du,du1,du2,k1,k2,k3,k4,
                    fsalfirst,fsallast,dT,nlsolver,tmp,atmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

struct Rosenbrock34ConstantCache{TF,UF,Tab,N,F} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
  nlsolver::N
  linsolve::F
end

function alg_cache(alg::Rodas3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  J,W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  linsolve = alg.linsolve(Val{:init},uf,u)
  nlsolver = SemiImplicitNLSolver(W,J,nothing,uf,nothing)
  Rosenbrock34ConstantCache(tf,uf,Rodas3Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits)),nlsolver,linsolve)
end

################################################################################

### ROS34PW methods

@ROS34PW(:cache)

################################################################################

### ROS4 methods

@Rosenbrock4(:cache)
jac_cache(c::Rosenbrock4Cache) = (c.J,c.W)

###############################################################################

### Rodas methods

struct Rodas4ConstantCache{TF,UF,Tab,N,F} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
  nlsolver::N
  linsolve::F
end

@cache mutable struct Rodas4Cache{uType,rateType,uNoUnitsType,N,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
  u::uType
  uprev::uType
  dense1::rateType
  dense2::rateType
  du::rateType
  du1::rateType
  du2::rateType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  fsalfirst::rateType
  fsallast::rateType
  dT::rateType
  nlsolver::N
  tmp::rateType
  atmp::uNoUnitsType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

function alg_cache(alg::Rodas4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  dense1 = zero(rate_prototype)
  dense2 = zero(rate_prototype)
  du = zero(rate_prototype)
  du1 = zero(rate_prototype)
  du2 = zero(rate_prototype)
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k6 = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  fsallast = zero(rate_prototype)
  dT = zero(rate_prototype)
  J,W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  tmp = zero(rate_prototype)
  atmp = similar(u, uEltypeNoUnits)
  tab = Rodas4Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  nlsolver = SemiImplicitNLSolver(W,J,du1,uf,jac_config)
  Rodas4Cache(u,uprev,dense1,dense2,du,du1,du2,k1,k2,k3,k4,
                    k5,k6,
                    fsalfirst,fsallast,dT,nlsolver,tmp,atmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Rodas4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  J,W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  linsolve = alg.linsolve(Val{:init},uf,u)
  nlsolver = SemiImplicitNLSolver(W,J,nothing,uf,nothing)
  Rodas4ConstantCache(tf,uf,Rodas4Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits)),nlsolver,linsolve)
end

function alg_cache(alg::Rodas42,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  dense1 = zero(rate_prototype)
  dense2 = zero(rate_prototype)
  du = zero(rate_prototype)
  du1 = zero(rate_prototype)
  du2 = zero(rate_prototype)
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k6 = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  fsallast = zero(rate_prototype)
  dT = zero(rate_prototype)
  J,W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  tmp = zero(rate_prototype)
  atmp = similar(u, uEltypeNoUnits)
  tab = Rodas42Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  nlsolver = SemiImplicitNLSolver(W,J,du1,uf,jac_config)
  Rodas4Cache(u,uprev,dense1,dense2,du,du1,du2,k1,k2,k3,k4,
                    k5,k6,
                    fsalfirst,fsallast,dT,nlsolver,tmp,atmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Rodas42,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  J,W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  linsolve = alg.linsolve(Val{:init},uf,u)
  nlsolver = SemiImplicitNLSolver(W,J,nothing,uf,nothing)
  Rodas4ConstantCache(tf,uf,Rodas42Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits)),nlsolver,linsolve)
end

function alg_cache(alg::Rodas4P,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  dense1 = zero(rate_prototype)
  dense2 = zero(rate_prototype)
  du = zero(rate_prototype)
  du1 = zero(rate_prototype)
  du2 = zero(rate_prototype)
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k6 = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  fsallast = zero(rate_prototype)
  dT = zero(rate_prototype)
  J,W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  tmp = zero(rate_prototype)
  atmp = similar(u, uEltypeNoUnits)
  tab = Rodas4PTableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  nlsolver = SemiImplicitNLSolver(W,J,du1,uf,jac_config)
  Rodas4Cache(u,uprev,dense1,dense2,du,du1,du2,k1,k2,k3,k4,
                    k5,k6,
                    fsalfirst,fsallast,dT,nlsolver,tmp,atmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Rodas4P,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  J,W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  linsolve = alg.linsolve(Val{:init},uf,u)
  nlsolver = SemiImplicitNLSolver(W,J,nothing,uf,nothing)
  Rodas4ConstantCache(tf,uf,Rodas4PTableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits)),nlsolver,linsolve)
end

################################################################################

### Rosenbrock5

struct Rosenbrock5ConstantCache{TF,UF,Tab,N,F} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
  nlsolver::N
  linsolve::F
end

@cache mutable struct Rosenbrock5Cache{uType,rateType,uNoUnitsType,N,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
  u::uType
  uprev::uType
  dense1::rateType
  dense2::rateType
  du::rateType
  du1::rateType
  du2::rateType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  fsalfirst::rateType
  fsallast::rateType
  dT::rateType
  nlsolver::N
  tmp::rateType
  atmp::uNoUnitsType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

function alg_cache(alg::Rodas5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  dense1 = zero(rate_prototype)
  dense2 = zero(rate_prototype)
  du = zero(rate_prototype)
  du1 = zero(rate_prototype)
  du2 = zero(rate_prototype)
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k6 = zero(rate_prototype)
  k7 = zero(rate_prototype)
  k8 = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  fsallast = zero(rate_prototype)
  dT = zero(rate_prototype)
  J,W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  tmp = zero(rate_prototype)
  atmp = similar(u, uEltypeNoUnits)
  tab = Rodas5Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  nlsolver = SemiImplicitNLSolver(W,J,du1,uf,jac_config)
  Rosenbrock5Cache(u,uprev,dense1,dense2,du,du1,du2,k1,k2,k3,k4,
                    k5,k6,k7,k8,
                    fsalfirst,fsallast,dT,nlsolver,tmp,atmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Rodas5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  J,W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  linsolve = alg.linsolve(Val{:init},uf,u)
  nlsolver = SemiImplicitNLSolver(W,J,nothing,uf,nothing)
  Rosenbrock5ConstantCache(tf,uf,Rodas5Tableau(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits)),nlsolver,linsolve)
end

################################################################################

### RosenbrockW6S4O

@RosenbrockW6S4OS(:cache)
