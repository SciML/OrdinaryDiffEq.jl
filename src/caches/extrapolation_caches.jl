@cache mutable struct AitkenNevilleCache{uType,rateType,arrayType,dtType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
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
  u_tmps::Array{uType,1}
  k_tmps::Array{rateType,1}
end

@cache mutable struct AitkenNevilleConstantCache{dtType,arrayType} <: OrdinaryDiffEqConstantCache
  dtpropose::dtType
  T::arrayType
  cur_order::Int
  work::dtType
  A::Int
  step_no::Int
end

function alg_cache(alg::AitkenNeville,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  utilde = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  cur_order = max(alg.init_order, alg.min_order)
  dtpropose = zero(dt)
  T = Array{typeof(u),2}(undef, alg.max_order, alg.max_order)
  # Array of arrays of length equal to number of threads to store intermediate
  # values of u and k. [Thread Safety]
  u_tmps = Array{typeof(u),1}(undef, Threads.nthreads())
  k_tmps = Array{typeof(k),1}(undef, Threads.nthreads())
  # Initialize each element of u_tmps and k_tmps to different instance of
  # zeros array similar to u and k respectively
  for i=1:Threads.nthreads()
      u_tmps[i] = zero(u)
      k_tmps[i] = zero(rate_prototype)
  end
  # Initialize lower triangle of T to different instance of zeros array similar to u
  for i=1:alg.max_order
    for j=1:i
      T[i,j] = zero(u)
    end
  end
  work = zero(dt)
  A = one(Int)
  atmp = similar(u,uEltypeNoUnits)
  step_no = zero(Int)
  AitkenNevilleCache(u,uprev,tmp,k,utilde,atmp,fsalfirst,dtpropose,T,cur_order,work,A,step_no,u_tmps,k_tmps)
end

function alg_cache(alg::AitkenNeville,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  dtpropose = zero(dt)
  cur_order = max(alg.init_order, alg.min_order)
  T = Array{typeof(u),2}(undef, alg.max_order, alg.max_order)
  @.. T = u
  work = zero(dt)
  A = one(Int)
  step_no = zero(Int)
  AitkenNevilleConstantCache(dtpropose,T,cur_order,work,A,step_no)
end


struct extrapolation_coefficients
  # This structure is used by the caches of the algorithms
  # ExtrapolationMidpointDeuflhard() and  ExtrapolationMidpointHairerWanner().
  # It contains the constant coefficients used to extrapolate the internal discretisations
  # in their perfom_step! function and some additional constant data.

  subdividing_sequence::Array{BigInt,1}  # subdividing_sequence[n] is used for the (n -1)th internal discretisation

  # Weights and Scaling factors for extrapolation operators
  extrapolation_weights::Array{Rational{BigInt},2}
  extrapolation_scalars::Array{Rational{BigInt},1}

  # Weights and scaling factors for internal extrapolation operators (used for error estimate)
  extrapolation_weights_2::Array{Rational{BigInt},2}
  extrapolation_scalars_2::Array{Rational{BigInt},1}
end

function create_extrapolation_coefficients(alg::algType) where {algType <: Union{ExtrapolationMidpointDeuflhard, ExtrapolationMidpointHairerWanner}}
  # Compute and return extrapolation_coefficients

  @unpack n_min, n_init, n_max, sequence = alg

  # Initialize subdividing_sequence:
  if sequence == :harmonic
      subdividing_sequence = [BigInt(n+1) for n = 0:n_max]
  elseif sequence == :romberg
      subdividing_sequence = [BigInt(2)^n for n = 0:n_max]
  else # sequence == :bulirsch
      subdividing_sequence = [n==0 ? BigInt(1) : (isodd(n) ? BigInt(2)^Int64(n/2 + 0.5) : 3BigInt(2^Int64(n/2 - 1))) for n = 0:n_max]
  end


  # Compute nodes corresponding to subdividing_sequence
  nodes = BigInt(1) .// subdividing_sequence .^ 2

  # Compute barycentric weights for internal extrapolation operators
  extrapolation_weights_2 = zeros(Rational{BigInt}, n_max, n_max)
  extrapolation_weights_2[1,:] = ones(Rational{BigInt}, 1, n_max)
  for n = 2:n_max
      distance = nodes[2:n] .- nodes[n+1]
      extrapolation_weights_2[1:(n-1), n] = extrapolation_weights_2[1:n-1, n-1] .// distance
      extrapolation_weights_2[n, n] = 1 // prod(-distance)
  end

  # Compute barycentric weights for extrapolation operators
  extrapolation_weights = zeros(Rational{BigInt}, n_max+1, n_max+1)
  for n = 1:n_max
      extrapolation_weights[n+1, (n+1) : (n_max+1)] = extrapolation_weights_2[n, n:n_max] // (nodes[n+1] - nodes[1])
      extrapolation_weights[1, n] = 1 // prod(nodes[1] .- nodes[2:n])
  end
  extrapolation_weights[1, n_max+1] = 1 // prod(nodes[1] .- nodes[2:n_max+1])

  # Rescale barycentric weights to obtain weights of 1. Barycentric Formula
  for m = 1:(n_max+1)
      extrapolation_weights[1:m, m] = - extrapolation_weights[1:m, m] .// nodes[1:m]
      if 2 <= m
          extrapolation_weights_2[1:m-1, m-1] = -extrapolation_weights_2[1:m-1, m-1] .// nodes[2:m]
      end
  end

  # Compute scaling factors for internal extrapolation operators
  extrapolation_scalars_2 = ones(Rational{BigInt}, n_max)
  extrapolation_scalars_2[1] = -nodes[2]
  for n = 1:(n_max-1)
      extrapolation_scalars_2[n+1] = -extrapolation_scalars_2[n] * nodes[n+2]
  end

  # Compute scaling factors for extrapolation operators
  extrapolation_scalars = -nodes[1] * [BigInt(1); extrapolation_scalars_2]

  # Initialize and return extrapolation_coefficients
  extrapolation_coefficients(subdividing_sequence,
      extrapolation_weights, extrapolation_scalars,
      extrapolation_weights_2, extrapolation_scalars_2)
end

@cache mutable struct ExtrapolationMidpointDeuflhardConstantCache{QType, extrapolation_coefficients} <: OrdinaryDiffEqConstantCache
  # Values that are mutated
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n + alg.n_min - 1)
  n_curr::Int64 # Storage for the current extrapolation order
  n_old::Int64 # Storage for the extrapolation order n_curr before perfom_step! changes the latter

  # Constant values
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int64} # stage_number[n] contains information for extrapolation order (n + alg.n_min - 1)
end

function alg_cache(alg::ExtrapolationMidpointDeuflhard,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
    # Initialize cache's members
    QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl

    Q = fill(zero(QType),alg.n_max - alg.n_min + 1)
    n_curr = alg.n_init
    n_old = alg.n_init

    coefficients = create_extrapolation_coefficients(alg)
    stage_number = [2sum(Int64.(coefficients.subdividing_sequence[1:n+1])) - n for n = alg.n_min:alg.n_max]

    # Initialize cache
    ExtrapolationMidpointDeuflhardConstantCache(Q, n_curr, n_old, coefficients, stage_number)
end



@cache mutable struct ExtrapolationMidpointDeuflhardCache{uType,uNoUnitsType,rateType,QType,extrapolation_coefficients} <: OrdinaryDiffEqMutableCache
  # Values that are mutated
  utilde::uType
  u_temp1::uType
  u_temp2::uType
  u_temp3::Array{uType,1}
  u_temp4::Array{uType,1}
  tmp::uType # for get_tmp_cache()
  T::Array{uType,1}  # Storage for the internal discretisations obtained by the explicit midpoint rule
  res::uNoUnitsType # Storage for the scaled residual of u and utilde

  fsalfirst::rateType
  k::rateType
  k_tmps::Array{rateType,1}

  # Constant values
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n + alg.n_min - 1)
  n_curr::Int64 # Storage for the current extrapolation order
  n_old::Int64 # Storage for the extrapolation order n_curr before perfom_step! changes the latter
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int64} # Stage_number[n] contains information for extrapolation order (n + alg.n_min - 1)
end

function alg_cache(alg::ExtrapolationMidpointDeuflhard,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  # Initialize cache's members
  utilde = zero(u)
  u_temp1 = zero(u)
  u_temp2 = zero(u)
  u_temp3 = Array{typeof(u),1}(undef, Threads.nthreads())
  u_temp4 = Array{typeof(u),1}(undef, Threads.nthreads())
  
  for i=1:Threads.nthreads()
      u_temp3[i] = zero(u)
      u_temp4[i] = zero(u)
  end

  tmp = zero(u)
  T = fill(zero(u), alg.n_max + 1)
  res = uEltypeNoUnits.(zero(u))

  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  k_tmps = Array{typeof(k),1}(undef, Threads.nthreads())
  for i=1:Threads.nthreads()
      k_tmps[i] = zero(rate_prototype)
  end

  cc =  alg_cache(alg::ExtrapolationMidpointDeuflhard,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,Val{false})
  # Initialize cache
  ExtrapolationMidpointDeuflhardCache(utilde, u_temp1, u_temp2, u_temp3, u_temp4, tmp, T, res, fsalfirst, k, k_tmps, cc.Q, cc.n_curr, cc.n_old, cc.coefficients,cc.stage_number)
end

@cache mutable struct ExtrapolationMidpointHairerWannerConstantCache{QType,extrapolation_coefficients} <: OrdinaryDiffEqConstantCache
  # Values that are mutated
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
  n_curr::Int64 # Storage for the current extrapolation order
  n_old::Int64 # Storage for the extrapolation order n_curr before perfom_step! changes the latter

  # Constant values
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int64} # stage_number[n] contains information for extrapolation order (n - 1)
  sigma::Rational{Int64} # Parameter for order selection
end

function alg_cache(alg::ExtrapolationMidpointHairerWanner,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  # Initialize cache's members
  QType = tTypeNoUnits <: Integer ? typeof(qmin_default(alg)) : tTypeNoUnits # Cf. DiffEqBase.__init in solve.jl

  Q = fill(zero(QType),alg.n_max + 1)
  n_curr = alg.n_init
  n_old = alg.n_init

  coefficients = create_extrapolation_coefficients(alg)
  stage_number = [2sum(Int64.(coefficients.subdividing_sequence[1:n+1])) - n for n = 0:alg.n_max]
  sigma = 9//10

  # Initialize the constant cache
  ExtrapolationMidpointHairerWannerConstantCache(Q, n_curr, n_old, coefficients, stage_number, sigma)
end

@cache mutable struct ExtrapolationMidpointHairerWannerCache{uType,uNoUnitsType,rateType,QType,extrapolation_coefficients} <: OrdinaryDiffEqMutableCache
  # Values that are mutated
  utilde::uType
  u_temp1::uType
  u_temp2::uType
  tmp::uType # for get_tmp_cache()
  T::Array{uType,1}  # Storage for the internal discretisations obtained by the explicit midpoint rule
  res::uNoUnitsType # Storage for the scaled residual of u and utilde

  fsalfirst::rateType
  k::rateType

  # Constant values
  Q::Vector{QType} # Storage for stepsize scaling factors. Q[n] contains information for extrapolation order (n - 1)
  n_curr::Int64 # Storage for the current extrapolation order
  n_old::Int64 # Storage for the extrapolation order n_curr before perfom_step! changes the latter
  coefficients::extrapolation_coefficients
  stage_number::Vector{Int64} # stage_number[n] contains information for extrapolation order (n - 1)
  sigma::Rational{Int64} # Parameter for order selection
end


function alg_cache(alg::ExtrapolationMidpointHairerWanner,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  # Initialize cache's members
  utilde = zero(u)
  u_temp1 = zero(u)
  u_temp2 = zero(u)
  tmp = zero(u)
  T = fill(zero(u), alg.n_max + 1)
  res = uEltypeNoUnits.(zero(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)

  cc = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,Val{false})

  # Initialize the cache
  ExtrapolationMidpointHairerWannerCache(utilde, u_temp1, u_temp2, tmp, T, res, fsalfirst, k,
      cc.Q, cc.n_curr, cc.n_old, cc.coefficients, cc.stage_number, cc.sigma)
end
