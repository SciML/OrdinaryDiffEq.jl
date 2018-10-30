abstract type RosenbrockMutableCache <: OrdinaryDiffEqMutableCache end
################################################################################

# Shampine's Low-order Rosenbrocks

mutable struct Rosenbrock23Cache{uType,rateType,JType,WType,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
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
  J::JType
  W::WType
  tmp::rateType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

u_cache(c::Rosenbrock23Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock23Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock23Cache) = (c.J,c.W)

mutable struct Rosenbrock32Cache{uType,rateType,JType,WType,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
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
  J::JType
  W::WType
  tmp::rateType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

u_cache(c::Rosenbrock32Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock32Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock32Cache) = (c.J,c.W)
user_cache(cache::Union{Rosenbrock23Cache,Rosenbrock32Cache}) = (cache.u,cache.uprev,cache.jac_config.duals...,cache.grad_config.duals)

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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = Rosenbrock23ConstantCache(uEltypeNoUnits,identity,identity)
  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)

  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)

  Rosenbrock23Cache(u,uprev,k₁,k₂,k₃,du1,du2,f₁,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = Rosenbrock32ConstantCache(uEltypeNoUnits,identity,identity)

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rosenbrock32Cache(u,uprev,k₁,k₂,k₃,du1,du2,f₁,fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,linsolve,jac_config,grad_config)
end

struct Rosenbrock23ConstantCache{T,TF,UF} <: OrdinaryDiffEqConstantCache
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

function alg_cache(alg::Rosenbrock23,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rosenbrock23ConstantCache(uEltypeNoUnits,tf,uf)
end

struct Rosenbrock32ConstantCache{T,TF,UF} <: OrdinaryDiffEqConstantCache
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

function alg_cache(alg::Rosenbrock32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rosenbrock32ConstantCache(uEltypeNoUnits,tf,uf)
end

################################################################################

### 3rd order specialized Rosenbrocks

struct Rosenbrock33ConstantCache{TF,UF,Tab} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
end

mutable struct Rosenbrock33Cache{uType,rateType,JType,WType,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
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
  J::JType
  W::WType
  tmp::rateType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

u_cache(c::Rosenbrock33Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock33Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock33Cache) = (c.J,c.W)

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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = ROS3PConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rosenbrock33Cache(u,uprev,du,du1,du2,k1,k2,k3,k4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::ROS3P,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rosenbrock33ConstantCache(tf,uf,ROS3PConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
end

mutable struct Rosenbrock34Cache{uType,rateType,JType,WType,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
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
  J::JType
  W::WType
  tmp::rateType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

u_cache(c::Rosenbrock34Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock34Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock34Cache) = (c.J,c.W)

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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = Rodas3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rosenbrock34Cache(u,uprev,du,du1,du2,k1,k2,k3,k4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

struct Rosenbrock34ConstantCache{TF,UF,Tab} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
end

function alg_cache(alg::Rodas3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rosenbrock34ConstantCache(tf,uf,Rodas3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
end


################################################################################

struct Rosenbrock4ConstantCache{TF,UF,Tab} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
end

mutable struct Rosenbrock4Cache{uType,rateType,JType,WType,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
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
  J::JType
  W::WType
  tmp::rateType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

u_cache(c::Rosenbrock4Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock4Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock4Cache) = (c.J,c.W)

function alg_cache(alg::RosShamp4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = RosShamp4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rosenbrock4Cache(u,uprev,du,du1,du2,k1,k2,k3,k4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::RosShamp4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rosenbrock4ConstantCache(tf,uf,RosShamp4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
end

function alg_cache(alg::Veldd4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = Veldd4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rosenbrock4Cache(u,uprev,du,du1,du2,k1,k2,k3,k4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Veldd4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rosenbrock4ConstantCache(tf,uf,Veldd4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
end

function alg_cache(alg::Velds4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = Velds4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rosenbrock4Cache(u,uprev,du,du1,du2,k1,k2,k3,k4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Velds4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rosenbrock4ConstantCache(tf,uf,Velds4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
end

function alg_cache(alg::GRK4T,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = GRK4TConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rosenbrock4Cache(u,uprev,du,du1,du2,k1,k2,k3,k4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::GRK4T,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rosenbrock4ConstantCache(tf,uf,GRK4TConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
end

function alg_cache(alg::GRK4A,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = GRK4AConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rosenbrock4Cache(u,uprev,du,du1,du2,k1,k2,k3,k4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::GRK4A,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rosenbrock4ConstantCache(tf,uf,GRK4AConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
end

function alg_cache(alg::Ros4LStab,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = Ros4LStabConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rosenbrock4Cache(u,uprev,du,du1,du2,k1,k2,k3,k4,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Ros4LStab,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rosenbrock4ConstantCache(tf,uf,Ros4LStabConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
end

###############################################################################

### Rodas methods

struct Rodas4ConstantCache{TF,UF,Tab} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
end

mutable struct Rodas4Cache{uType,rateType,JType,WType,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
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
  J::JType
  W::WType
  tmp::rateType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

u_cache(c::Rodas4Cache) = (c.dT,c.tmp)
du_cache(c::Rodas4Cache) = (c.dense1,c.dense2,c.du,c.du1,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rodas4Cache) = (c.J,c.W)

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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = Rodas4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rodas4Cache(u,uprev,dense1,dense2,du,du1,du2,k1,k2,k3,k4,
                    k5,k6,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Rodas4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rodas4ConstantCache(tf,uf,Rodas4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = Rodas42ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rodas4Cache(u,uprev,dense1,dense2,du,du1,du2,k1,k2,k3,k4,
                    k5,k6,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Rodas42,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rodas4ConstantCache(tf,uf,Rodas42ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = Rodas4PConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rodas4Cache(u,uprev,dense1,dense2,du,du1,du2,k1,k2,k3,k4,
                    k5,k6,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Rodas4P,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rodas4ConstantCache(tf,uf,Rodas4PConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
end

################################################################################

### Rosenbrock5

struct Rosenbrock5ConstantCache{TF,UF,Tab} <: OrdinaryDiffEqConstantCache
  tf::TF
  uf::UF
  tab::Tab
end

mutable struct Rosenbrock5Cache{uType,rateType,JType,WType,TabType,TFType,UFType,F,JCType,GCType} <: RosenbrockMutableCache
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
  J::JType
  W::WType
  tmp::rateType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::rateType
  linsolve::F
  jac_config::JCType
  grad_config::GCType
end

u_cache(c::Rosenbrock5Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock5Cache) = (c.dense1,c.dense2,c.du,c.du1,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock5Cache) = (c.J,c.W)

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
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  tmp = zero(rate_prototype)
  tab = Rodas5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  tf = DiffEqDiffTools.TimeGradientWrapper(f,uprev,p)
  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve_tmp = zero(rate_prototype)
  linsolve = alg.linsolve(Val{:init},uf,u)
  grad_config = build_grad_config(alg,f,tf,du1,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  Rosenbrock5Cache(u,uprev,dense1,dense2,du,du1,du2,k1,k2,k3,k4,
                    k5,k6,k7,k8,
                    fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,
                    linsolve,jac_config,grad_config)
end

function alg_cache(alg::Rodas5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tf = DiffEqDiffTools.TimeDerivativeWrapper(f,u,p)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  Rosenbrock5ConstantCache(tf,uf,Rodas5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits)))
end
