mutable struct RadauIIA5ConstantCache{F,Tab,Tol,Dt,U,JType} <: OrdinaryDiffEqConstantCache
  uf::F
  tab::Tab
  κ::Tol
  ηold::Tol
  iter::Int
  cont1::U
  cont2::U
  cont3::U
  dtprev::Dt
  W_γdt::Dt
  status::DiffEqBase.NLStatus
  J::JType
end

function alg_cache(alg::RadauIIA5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  uf  = UDerivativeWrapper(f, t, p)
  uToltype = real(uBottomEltypeNoUnits)
  tab = RadauIIA5Tableau(uToltype, real(tTypeNoUnits))

  κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1//100)
  J = false .* _vec(rate_prototype) .* _vec(rate_prototype)'

  RadauIIA5ConstantCache(uf, tab, κ, one(uToltype), 10000, u, u, u, dt, dt, DiffEqBase.Convergence, J)
end

mutable struct RadauIIA5Cache{uType,cuType,uNoUnitsType,rateType,JType,W1Type,W2Type,UF,JC,F1,F2,Tab,Tol,Dt,rTol,aTol} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  z1::uType
  z2::uType
  z3::uType
  w1::uType
  w2::uType
  w3::uType
  dw1::uType
  dw23::cuType
  cont1::uType
  cont2::uType
  cont3::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  k2::rateType
  k3::rateType
  fw1::rateType
  fw2::rateType
  fw3::rateType
  J::JType
  W1::W1Type
  W2::W2Type # complex
  uf::UF
  tab::Tab
  κ::Tol
  ηold::Tol
  iter::Int
  tmp::uType
  atmp::uNoUnitsType
  jac_config::JC
  linsolve1::F1
  linsolve2::F2
  rtol::rTol
  atol::aTol
  dtprev::Dt
  W_γdt::Dt
  status::DiffEqBase.NLStatus
end

function alg_cache(alg::RadauIIA5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  uf = UJacobianWrapper(f, t, p)
  uToltype = real(uBottomEltypeNoUnits)
  tab = RadauIIA5Tableau(uToltype, real(tTypeNoUnits))

  κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1//100)

  z1 = similar(u); z2 = similar(u); z3 = similar(u)
  w1 = similar(u); w2 = similar(u); w3 = similar(u)
  dw1 = similar(u); dw23 = similar(u, Complex{eltype(u)})
  cont1 = similar(u); cont2 = similar(u); cont3 = similar(u)

  fsalfirst = similar(rate_prototype)
  k = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype)
  fw1 = similar(rate_prototype); fw2 = similar(rate_prototype); fw3 = similar(rate_prototype)

  J = false .* vec(rate_prototype) .* vec(rate_prototype)'
  W1 = similar(J); W2 = similar(J, Complex{eltype(J)})

  du1 = similar(rate_prototype)

  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dw1)
  linsolve1 = alg.linsolve(Val{:init}, uf, u)
  linsolve2 = alg.linsolve(Val{:init}, uf, u)
  rtol = reltol isa Number ? reltol : similar(reltol)
  atol = reltol isa Number ? reltol : similar(reltol)

  RadauIIA5Cache(u, uprev,
                 z1, z2, z3, w1, w2, w3,
                 dw1, dw23, cont1, cont2, cont3,
                 du1, fsalfirst, k, k2, k3, fw1, fw2, fw3,
                 J, W1, W2,
                 uf, tab, κ, one(uToltype), 10000,
                 tmp, atmp, jac_config, linsolve1, linsolve2, rtol, atol, dt, dt, DiffEqBase.Convergence)
end
