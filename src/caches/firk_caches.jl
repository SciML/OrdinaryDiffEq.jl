mutable struct RadauIIA5ConstantCache{F,Tab,Tol,Dt,U} <: OrdinaryDiffEqConstantCache
  uf::F
  tab::Tab
  κ::Tol
  tol::Tol
  ηold::Tol
  nl_iters::Int
  dtprev::Dt
  cont1::U
  cont2::U
  cont3::U
end

function alg_cache(alg::RadauIIA5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})

  uf  = DiffEqDiffTools.UDerivativeWrapper(f, t, p)
  uToltype = real(uBottomEltypeNoUnits)
  tab = RadauIIA5Tableau(uToltype, real(tTypeNoUnits))

  κ = alg.κ !== nothing ? uToltype(alg.κ) : uToltype(1//100)
  tol = alg.tol !== nothing ? uToltype(alg.tol) : uToltype(min(0.03,first(reltol)^(0.5)))

  RadauIIA5ConstantCache(uf, tab, κ, tol, zero(tol), 10000, dt, u, u, u)
end
