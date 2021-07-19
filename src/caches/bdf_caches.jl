@cache mutable struct ABDF2ConstantCache{N,dtType,rate_prototype} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  eulercache::ImplicitEulerConstantCache
  dtₙ₋₁::dtType
  fsalfirstprev::rate_prototype
end

function alg_cache(alg::ABDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  γ, c = 2//3, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  eulercache = ImplicitEulerConstantCache(nlsolver)

  dtₙ₋₁ = one(dt)
  fsalfirstprev = rate_prototype

  ABDF2ConstantCache(nlsolver, eulercache, dtₙ₋₁, fsalfirstprev)
end

@cache mutable struct ABDF2Cache{uType,rateType,uNoUnitsType,N,dtType} <: OrdinaryDiffEqMutableCache
  uₙ::uType
  uₙ₋₁::uType
  uₙ₋₂::uType
  fsalfirst::rateType
  fsalfirstprev::rateType
  zₙ₋₁::uType
  atmp::uNoUnitsType
  nlsolver::N
  eulercache::ImplicitEulerCache
  dtₙ₋₁::dtType
end

function alg_cache(alg::ABDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  γ, c = 2//3, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  fsalfirstprev = zero(rate_prototype)
  atmp = similar(u,uEltypeNoUnits)

  eulercache = ImplicitEulerCache(u,uprev,uprev2,fsalfirst,atmp,nlsolver)

  dtₙ₋₁ = one(dt)
  zₙ₋₁ = zero(u)

  ABDF2Cache(u,uprev,uprev2,fsalfirst,fsalfirstprev,zₙ₋₁,atmp,
              nlsolver,eulercache,dtₙ₋₁)
end

# SBDF

@cache mutable struct SBDFConstantCache{rateType,N,uType} <: OrdinaryDiffEqConstantCache
  cnt::Int
  k2::rateType
  nlsolver::N
  uprev2::uType
  uprev4::uType
  uprev3::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du₁::rateType
  du₂::rateType
end

@cache mutable struct SBDFCache{uType,rateType,N} <: OrdinaryDiffEqMutableCache
  cnt::Int
  u::uType
  uprev::uType
  fsalfirst::rateType
  nlsolver::N
  uprev2::uType
  uprev3::uType
  uprev4::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du₁::rateType
  du₂::rateType
end

function alg_cache(alg::SBDF,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  γ, c = 1//1, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))

  k2 = rate_prototype
  k₁ = rate_prototype; k₂ = rate_prototype; k₃ = rate_prototype
  du₁ = rate_prototype; du₂ = rate_prototype

  uprev2 = u; uprev3 = u; uprev4 = u

  SBDFConstantCache(1,k2,nlsolver,uprev2,uprev3,uprev4,k₁,k₂,k₃,du₁,du₂)
end

function alg_cache(alg::SBDF,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  γ, c = 1//1, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  order = alg.order

  k₁ = zero(rate_prototype)
  k₂ = order >= 3 ? zero(rate_prototype) : k₁
  k₃ = order == 4 ? zero(rate_prototype) : k₁
  du₁ = zero(rate_prototype)
  du₂ = zero(rate_prototype)

  uprev2 = zero(u)
  uprev3 = order >= 3 ? zero(u) : uprev2
  uprev4 = order == 4 ? zero(u) : uprev2

  SBDFCache(1,u,uprev,fsalfirst,nlsolver,uprev2,uprev3,uprev4,k₁,k₂,k₃,du₁,du₂)
end

# QNDF1

@cache mutable struct QNDF1ConstantCache{N,coefType,coefType1,coefType2,dtType,uType} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  D::coefType1
  D2::coefType2
  R::coefType
  U::coefType
  uprev2::uType
  dtₙ₋₁::dtType
end

@cache mutable struct QNDF1Cache{uType,rateType,coefType,coefType1,coefType2,uNoUnitsType,N,dtType} <: OrdinaryDiffEqMutableCache
  uprev2::uType
  fsalfirst::rateType
  D::coefType1
  D2::coefType2
  R::coefType
  U::coefType
  atmp::uNoUnitsType
  utilde::uType
  nlsolver::N
  dtₙ₋₁::dtType
end

function alg_cache(alg::QNDF1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  γ, c = zero(inv((1-alg.kappa))), 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))

  uprev2 = u
  dtₙ₋₁ = zero(t)

  D = fill(zero(u), 1, 1)
  D2 = fill(zero(u), 1, 2)
  R = fill(zero(t), 1, 1)
  U = fill(zero(t), 1, 1)

  U!(1,U)

  QNDF1ConstantCache(nlsolver,D,D2,R,U,uprev2,dtₙ₋₁)
end

function alg_cache(alg::QNDF1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  γ, c = zero(inv((1-alg.kappa))), 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  D = Array{typeof(u)}(undef, 1, 1)
  D2 = Array{typeof(u)}(undef, 1, 2)

  R = fill(zero(t), 1, 1)
  U = fill(zero(t), 1, 1)

  D[1] = zero(u)
  D2[1] = zero(u); D2[2] = zero(u)

  U!(1,U)

  atmp = similar(u,uEltypeNoUnits)
  utilde = zero(u)
  uprev2 = zero(u)
  dtₙ₋₁ = zero(dt)

  QNDF1Cache(uprev2,fsalfirst,D,D2,R,U,atmp,utilde,nlsolver,dtₙ₋₁)
end

# QNDF2

@cache mutable struct QNDF2ConstantCache{N,coefType,coefType1,coefType2,uType,dtType} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  D::coefType1
  D2::coefType2
  R::coefType
  U::coefType
  uprev2::uType
  uprev3::uType
  dtₙ₋₁::dtType
  dtₙ₋₂::dtType
end

@cache mutable struct QNDF2Cache{uType,rateType,coefType,coefType1,coefType2,uNoUnitsType,N,dtType} <: OrdinaryDiffEqMutableCache
  uprev2::uType
  uprev3::uType
  fsalfirst::rateType
  D::coefType1
  D2::coefType2
  R::coefType
  U::coefType
  atmp::uNoUnitsType
  utilde::uType
  nlsolver::N
  dtₙ₋₁::dtType
  dtₙ₋₂::dtType
end

function alg_cache(alg::QNDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  γ, c = zero(inv((1-alg.kappa))), 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))

  uprev2 = u
  uprev3 = u
  dtₙ₋₁ = zero(t)
  dtₙ₋₂ = zero(t)

  D = fill(zero(u), 1, 2)
  D2 = fill(zero(u), 1, 3)
  R = fill(zero(t), 2, 2)
  U = fill(zero(t), 2, 2)

  U!(2,U)

  QNDF2ConstantCache(nlsolver,D,D2,R,U,uprev2,uprev3,dtₙ₋₁,dtₙ₋₂)
end

function alg_cache(alg::QNDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  γ, c = zero(inv((1-alg.kappa))), 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  D = Array{typeof(u)}(undef, 1, 2)
  D2 = Array{typeof(u)}(undef, 1, 3)
  R = fill(zero(t), 2, 2)
  U = fill(zero(t), 2, 2)

  D[1] = zero(u); D[2] = zero(u)
  D2[1] = zero(u);  D2[2] = zero(u); D2[3] = zero(u)

  U!(2,U)

  atmp = similar(u,uEltypeNoUnits)
  utilde = zero(u)
  uprev2 = zero(u)
  uprev3 = zero(u)
  dtₙ₋₁ = zero(dt)
  dtₙ₋₂ = zero(dt)

  QNDF2Cache(uprev2,uprev3,fsalfirst,D,D2,R,U,atmp,utilde,nlsolver,dtₙ₋₁,dtₙ₋₂)
end

@cache mutable struct QNDFConstantCache{MO,N,coefType,UType,dtType,EEstType,gammaType} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  U::UType
  D::coefType
  prevD::coefType
  prevorder::Int
  order::Int
  max_order::Val{MO}
  dtprev::dtType
  nconsteps::Int ##Successful Consecutive Step with the same step size
  consfailcnt::Int #Consecutive failed steps count
  EEst1::EEstType #Error Estimator for k-1 order
  EEst2::EEstType #Error Estimator for k+1 order
  γₖ::gammaType
end

function alg_cache(alg::QNDF{MO},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where MO
  max_order = MO
  γ, c = one(eltype(alg.kappa)), 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  dtprev = one(dt)
  D = Matrix{uEltypeNoUnits}(undef, length(u), max_order+2)
  fill!(D, zero(uEltypeNoUnits))
  prevD = similar(D)
  fill!(prevD, zero(uEltypeNoUnits))
  EEst1 = tTypeNoUnits(1)
  EEst2 = tTypeNoUnits(1)

  U = zero(MMatrix{max_order,max_order,tTypeNoUnits})
  for r = 1:max_order
    U[1,r] = -r
    for j = 2:max_order
      U[j,r] = U[j-1,r] * ((j-1) - r)/j
    end
  end
  U = SArray(U)

  γₖ = SVector(ntuple(k->sum(tTypeNoUnits(1//j) for j in 1:k), Val(max_order)))

  QNDFConstantCache(nlsolver, U, D, prevD, 1, 1, Val(max_order), dtprev, 0, 0, EEst1, EEst2, γₖ)
end

@cache mutable struct QNDFCache{MO,UType,RUType,rateType,N,coefType,dtType,EEstType,gammaType,uType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  fsalfirst::rateType
  dd::uType
  utilde::uType
  utildem1::uType
  utildep1::uType
  ϕ::uType
  u₀::uType
  nlsolver::N
  U::UType
  RU::RUType
  D::coefType
  Dtmp::coefType
  tmp2::uType
  prevD::coefType
  order::Int
  prevorder::Int
  max_order::Val{MO}
  dtprev::dtType
  nconsteps::Int ##Successful consecutive step with the same step size
  consfailcnt::Int #Consecutive failed steps count
  EEst1::EEstType #Error Estimator for k-1 order
  EEst2::EEstType #Error Estimator for k+1 order
  γₖ::gammaType
  atmp::uNoUnitsType
  atmpm1::uNoUnitsType
  atmpp1::uNoUnitsType
end

function alg_cache(alg::QNDF{MO},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where MO
  max_order = MO
  γ, c = one(eltype(alg.kappa)), 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)
  dd = zero(u)
  utilde = zero(u)
  utildem1 = zero(u)
  utildep1 = zero(u)
  ϕ = zero(u)
  u₀ = zero(u)
  dtprev = one(dt)
  D = similar(u, uEltypeNoUnits, length(u), max_order + 2)
  fill!(D, zero(uEltypeNoUnits))
  Dtmp = similar(D)
  fill!(Dtmp, zero(uEltypeNoUnits))
  prevD = zero(similar(D))
  atmp = zero(similar(u, uEltypeNoUnits))
  atmpm1 = zero(similar(u, uEltypeNoUnits))
  atmpp1 = zero(similar(u, uEltypeNoUnits))
  tmp2 = zero(u)
  EEst1 = tTypeNoUnits(1)
  EEst2 = tTypeNoUnits(1)

  U = zero(MMatrix{max_order,max_order,tTypeNoUnits})
  for r = 1:max_order
    U[1,r] = -r
    for j = 2:max_order
      U[j,r] = U[j-1,r] * ((j-1) - r)/j
    end
  end
  U = SArray(U)

  RU = Matrix(U)
  γₖ = SVector(ntuple(k->sum(tTypeNoUnits(1//j) for j in 1:k), Val(max_order)))

  QNDFCache(fsalfirst, dd, utilde, utildem1, utildep1, ϕ, u₀, nlsolver, U, RU, D, Dtmp, tmp2, prevD, 1, 1, Val(max_order), dtprev, 0, 0, EEst1, EEst2, γₖ, atmp, atmpm1, atmpp1)
end


@cache mutable struct MEBDF2Cache{uType,rateType,uNoUnitsType,N} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  tmp2::uType
  atmp::uNoUnitsType
  nlsolver::N
end

function alg_cache(alg::MEBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  γ, c = 1, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); tmp2 = zero(u)
  atmp = similar(u,uEltypeNoUnits)

  MEBDF2Cache(u,uprev,uprev2,fsalfirst,z₁,z₂,tmp2,atmp,nlsolver)
end

mutable struct MEBDF2ConstantCache{N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::MEBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  γ, c = 1, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  MEBDF2ConstantCache(nlsolver)
end

@cache mutable struct FBDFConstantCache{MO,N,tsType,tType,uType,uuType,coeffType,EEstType,rType,wType} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  ts::tsType
  ts_tmp::tsType
  t_old::tType
  u_history::uuType
  order::Int
  prev_order::Int
  u_corrector::uType
  bdf_coeffs::coeffType
  max_order::Val{MO}
  nconsteps::Int
  consfailcnt::Int
  terkm2::EEstType
  terkm1::EEstType
  terk::EEstType
  terkp1::EEstType
  r::rType
  weights::wType
  nonevesuccsteps::Int
end

function alg_cache(alg::FBDF{MO},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where MO
  γ, c = 1.0, 1.0
  max_order = MO
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  bdf_coeffs = SA[1 -1 0 0 0 0 ;
                  3//2 -2 1//2 0 0 0 ;
                  11//6 -3 3//2 -1//3  0 0 ;
                  25//12 -4 3 -4//3 1//4 0 ;
                  137//60 -5 5 -10//3 5//4 -1//5]
  ts = zero(Vector{typeof(t)}(undef,max_order+2)) #ts is the successful past points, it will be updated after successful step
  ts_tmp = similar(ts)

  u_history = zero(Matrix{eltype(u)}(undef,length(u),max_order+2))
  order = 1
  prev_order = 1
  u_corrector = similar(u_history)
  fill!(u_corrector,zero(eltype(u)))
  fill!(u_history,zero(eltype(u_history)))
  terkm2 = tTypeNoUnits(1)
  terkm1= tTypeNoUnits(1)
  terk= tTypeNoUnits(1)
  terkp1 = tTypeNoUnits(1)
  r = zero(Vector{typeof(t)}(undef,max_order+2)) 
  weights = zero(Vector{typeof(t)}(undef,max_order+2))
  weights[1] = 1
  nconsteps = 0
  consfailcnt = 0
  t_old = zero(t)
  nonevesuccsteps = 0
  
  FBDFConstantCache(nlsolver,ts,ts_tmp,t_old,u_history,order,prev_order,u_corrector,bdf_coeffs,Val(5),nconsteps,consfailcnt,terkm2,terkm1,terk,terkp1,r,weights,nonevesuccsteps)
end

@cache mutable struct FBDFCache{MO,N,rateType,uNoUnitsType,tsType,tType,uType,uuType,coeffType,EEstType,rType,wType} <: OrdinaryDiffEqMutableCache
  fsalfirst::rateType
  nlsolver::N
  ts::tsType
  ts_tmp::tsType
  t_old::tType
  u_history::uuType
  order::Int
  prev_order::Int
  u_corrector::uuType
  u₀::uType
  bdf_coeffs::coeffType
  max_order::Val{MO}
  nconsteps::Int
  consfailcnt::Int
  tmp::uType
  atmp::uNoUnitsType
  terkm2::EEstType
  terkm1::EEstType
  terk::EEstType
  terkp1::EEstType
  terk_tmp::uType
  terkp1_tmp::uType
  r::rType
  weights::wType
  equi_ts::tsType
  nonevesuccsteps::Int
end

function alg_cache(alg::FBDF{MO},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where MO
  γ, c = 1.0, 1.0
  fsalfirst = zero(rate_prototype)
  max_order = MO
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  bdf_coeffs = SA[1 -1 0 0 0 0 ;
                  3//2 -2 1//2 0 0 0 ;
                  11//6 -3 3//2 -1//3  0 0 ;
                  25//12 -4 3 -4//3 1//4 0 ;
                  137//60 -5 5 -10//3 5//4 -1//5]
  ts = zero(Vector{typeof(t)}(undef,max_order+2)) #ts is the successful past points, it will be updated after successful step
  u_history = zero(Matrix{eltype(u)}(undef,length(u),max_order+2))
  order = 1
  prev_order = 1
  u_corrector = similar(u_history)
  fill!(u_corrector,zero(eltype(u)))
  fill!(u_history,zero(eltype(u_history)))
  terkm2 = tTypeNoUnits(1)
  terkm1= tTypeNoUnits(1)
  terk= tTypeNoUnits(1)
  terkp1 = tTypeNoUnits(1)
  terk_tmp = similar(u)
  terkp1_tmp = similar(u)
  r = zero(Vector{typeof(t)}(undef,max_order+2)) 
  weights = zero(Vector{typeof(t)}(undef,max_order+2))
  weights[1] = 1
  nconsteps = 0
  consfailcnt = 0
  t_old = zero(t)
  atmp = similar(u, uEltypeNoUnits)
  fill!(atmp,zero(uEltypeNoUnits))
  u₀ = similar(u)
  equi_ts = similar(ts)
  tmp = similar(u)
  ts_tmp = similar(ts)
  nonevesuccsteps = 0

  FBDFCache(fsalfirst,nlsolver,ts,ts_tmp,t_old,u_history,order,prev_order,u_corrector,u₀,bdf_coeffs,Val(5),nconsteps,consfailcnt,tmp,atmp,terkm2,terkm1,terk,terkp1,terk_tmp,terkp1_tmp,r,weights,equi_ts,nonevesuccsteps)
end
