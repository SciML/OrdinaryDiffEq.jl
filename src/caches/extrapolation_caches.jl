@cache mutable struct RichardsonEulerCache{uType,rateType,arrayType,dtType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  utilde::uType
  atmp::uNoUnitsType
  fsalfirst::rateType
  dtpropose::dtType
  T::arrayType
  cur_order::Int
  work::dtType
  A::Int
  step_no::Int
end

@cache mutable struct RichardsonEulerConstantCache{dtType,arrayType} <: OrdinaryDiffEqConstantCache
  dtpropose::dtType
  T::arrayType
  cur_order::Int
  work::dtType
  A::Int
  step_no::Int
end

function alg_cache(alg::RichardsonEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  utilde = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  cur_order = max(alg.init_order, alg.min_order)
  dtpropose = zero(dt)
  T = fill(zeros(eltype(u), size(u)), (alg.max_order, alg.max_order))
  work = zero(dt)
  A = one(Int)
  atmp = similar(u,uEltypeNoUnits)
  step_no = zero(Int)
  RichardsonEulerCache(u,uprev,tmp,k,utilde,atmp,fsalfirst,dtpropose,T,cur_order,work,A,step_no)
end

function alg_cache(alg::RichardsonEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  dtpropose = zero(dt)
  cur_order = max(alg.init_order, alg.min_order)
  T = fill(zero(eltype(u)), (alg.max_order, alg.max_order))
  work = zero(dt)
  A = one(Int)
  step_no = zero(Int)
  RichardsonEulerConstantCache(dtpropose,T,cur_order,work,A,step_no)
end

@cache mutable struct ExtrapolationMidpointDeuflhardConstantCache{dtType,QType} <: OrdinaryDiffEqConstantCache
  dtpropose::dtType
  Q::Vector{QType} # storage for scaling factors of stepsize
  current_extrapolation_order::Int
  subdividing_sequence::Array{BigInt,1}
  # weights and scaling factors for extrapolation operators:
  extrapolation_weights::Array{Rational{BigInt},2}
  extrapolation_scalars::Array{Rational{BigInt},1}
  # weights and scaling factors for internal extrapolation operators (used for error estimate):
  extrapolation_weights_2::Array{Rational{BigInt},2}
  extrapolation_scalars_2::Array{Rational{BigInt},1}
end

function alg_cache(alg::ExtrapolationMidpointDeuflhard,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # cf. DiffEqBase.__init in solve.jl
  N = alg.max_extrapolation_order # for readability

  dtpropose = zero(dt)
  Q = fill(zero(QType),N - alg.min_extrapolation_order + 1)
  current_extrapolation_order = alg.init_extrapolation_order

  # initialize subdividing_sequence:
  if alg.sequence_symbol == :harmonic
      subdividing_sequence = [BigInt(n+1) for n = 0:N]
  elseif alg.sequence_symbol == :romberg
      subdividing_sequence = [BigInt(2)^n for n = 0:N]
  else # sequence_symbol == :bulirsch
      subdividing_sequence = [n==0 ? BigInt(1) : (isodd(n) ? BigInt(2)^Int64(n/2+0.5) : 3*BigInt(2^Int64(n/2-1))) for n = 0:N]
  end

  # compute nodes corresponding to the subdividing sequence subdividing_sequence
  nodes = BigInt(1).// subdividing_sequence.^2

  # compute barycentric weights for internal extrapolation operators
  extrapolation_weights_2 = zeros(Rational{BigInt},N,N)
  extrapolation_weights_2[1,:] = ones(Rational{BigInt},1,N)
  for n = 2:N
      distance = nodes[2:n] .- nodes[n+1]
      extrapolation_weights_2[1:(n-1),n] = extrapolation_weights_2[1:n-1,n-1] .// distance
      extrapolation_weights_2[n,n] = 1 // prod(-distance)
  end

  # compute barycentric weights for extrapolation operators
  extrapolation_weights = zeros(Rational{BigInt},N+1,N+1)
  for n = 1:N
      extrapolation_weights[n+1,(n+1):(N+1)] = extrapolation_weights_2[n,n:N] // (nodes[n+1]-nodes[1])
      extrapolation_weights[1,n] = 1 // prod(nodes[1].-nodes[2:n])
  end
  extrapolation_weights[1,N+1] = 1 // prod(nodes[1].-nodes[2:N+1])

  #rescale barycentric weights to obtain weights of 1. barycentric formula
  for m = 1:(N+1)
      extrapolation_weights[1:m,m] = - extrapolation_weights[1:m,m] .// nodes[1:m]
      if 2 <= m
          extrapolation_weights_2[1:m-1,m-1] = - extrapolation_weights_2[1:m-1,m-1] .// nodes[2:m]
      end
  end

  # compute scaling factors for internal extrapolation operators
  extrapolation_scalars_2 = ones(Rational{BigInt},N)
  extrapolation_scalars_2[1] = -nodes[2]
  for n = 1:(N-1)
      extrapolation_scalars_2[n+1] = -extrapolation_scalars_2[n]*nodes[n+2]
  end

  # compute scaling factors for extrapolation operators
  extrapolation_scalars = -nodes[1]*[BigInt(1); extrapolation_scalars_2]

  # initialize the constant cache
  ExtrapolationMidpointDeuflhardConstantCache(dtpropose, Q, current_extrapolation_order,subdividing_sequence, extrapolation_weights, extrapolation_scalars,extrapolation_weights_2, extrapolation_scalars_2)
end

@cache mutable struct ExtrapolationMidpointDeuflhardCache{uType,uNoUnitsType,rateType,dtType,QType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  utilde::uType
  tmp::uType
  atmp::uNoUnitsType
  k::rateType
  fsalfirst::rateType
  proposed_extrapolation_order::Int
  constant_cache::ExtrapolationMidpointDeuflhardConstantCache
end

function alg_cache(alg::ExtrapolationMidpointDeuflhard,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  utilde = similar(u)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  proposed_extrapolation_order = alg.init_extrapolation_order # order of first step is set by user
  constant_cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,Val{false})
  ExtrapolationMidpointDeuflhardCache(u,uprev,utilde,tmp,atmp,k,fsalfirst,proposed_extrapolation_order,constant_cache)
end
