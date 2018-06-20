const BIAS1 = 6
const BIAS2 = 6
const BIAS3 = 10
const ADDON = 1e-6

# This function computes the integral, from -1 to 0, of a polynomial
# `P(x)` from the coefficients of `P` with an offset `k`.
function ∫₋₁⁰dx(a, deg, k)
  @inbounds begin
    int = zero(eltype(a))
    sign = 1
    for i in 1:deg
      int += flipsign(a[i]/(i+k), sign)
      sign = -sign
    end
    return int
  end
end

# `l` is the coefficients of the polynomial `Λ` that satisfies conditions
# Λ(0) = 1, Λ(-1) = 0, and Λ̇(-ξᵢ) = 0, where ξᵢ = (tₙ-tₙ₋₁)/dt.
# It is described in the paper "A Polyalgorithm for the Numerical Solution
# of Ordinary Differential Equations" by G. D. Byrne and A. C. Hindmarsh in
# the page 86.
# https://dl.acm.org/citation.cfm?id=355636

# More implementation details are in the
# https://github.com/JuliaDiffEq/DiffEqDevMaterials repository
function calc_coeff!(cache::T) where T
  isvode = ( T <: JVODECache || T <: JVODEConstantCache )
  @inbounds begin
    isconst = T <: OrdinaryDiffEqConstantCache
    isconst || (cache = cache.const_cache)
    isvarorder = nordsieck_change_order(cache, 1)
    @unpack m, l, tau = cache
    dtsum = dt = tau[1]
    order = cache.step
    if order == 1
      l[1] = l[2] = cache.c_LTE₋₁ = cache.c_𝒟 = 1
      cache.c_LTE = 1//2
      cache.c_LTE₊₁ = 1//12
      cache.c_conv = 1//10 / cache.c_LTE
      return nothing
    end
    m[1] = 1
    for i in 2:order+1
      m[i] = 0
    end
    # initialize ξ_inv
    ξ_inv = dt / dtsum
    # compute coefficients from the Newton polynomial
    # check the `JuliaDiffEq/DiffEqDevMaterials` repository for more details
    for j in 1:order-1
      if isvarorder && j == order-1
        M₋₁ = ∫₋₁⁰dx(m, order-1, 1)
        # It is the same with `tq[1]` in SUNDIALS cvode.c
        cache.c_LTE₋₁ = order * M₋₁ / m[order-1]
      end
      ξ_inv = dt / dtsum
      for i in j:-1:1
        m[i+1] = muladd(m[i], ξ_inv, m[i+1])
      end
      dtsum += tau[j+1]
    end
    ξ_inv = dt / dtsum

    M0 = ∫₋₁⁰dx(m, order, 0)
    M1 = ∫₋₁⁰dx(m, order, 1)
    M0_inv = inv(M0)
    l[1] = 1
    for i in 1:order
      l[i+1] = M0_inv * m[i] / i
    end
    # TODO: simplify LTE calculation
    # This is the error estimation coefficient for the current order `q`
    # ||Δ||⋅c_LTE yields the difference between a `q` degree interpolating
    # polynomial and a `q+1` degree interpolating polynomial at time `t`.
    # It is the same with `tq[2]` in SUNDIALS cvode.c
    cache.c_LTE = M1 * M0_inv * ξ_inv
    # It is the same with `tq[5]` in SUNDIALS cvode.c
    isvode && (cache.c_𝒟 = inv(ξ_inv) / l[order+1])
    if isvarorder
      for i in order-1:-1:1
        m[i+1] = muladd(ξ_inv, m[i], m[i+1])
      end
      M2 = ∫₋₁⁰dx(m, order, 1)
      # It is the same with `tq[3]` in SUNDIALS cvode.c
      cache.c_LTE₊₁ = M2 * M0_inv / (order+1)
    end # endif isvarorder
    # It is the same with `tq[4]` in SUNDIALS cvode.c
    cache.c_conv = 1//10 / cache.c_LTE
    return nothing
  end # end @inbounds
end

# Apply the Pascal linear operator
function perform_predict!(cache::T, rewind=false) where T
  @inbounds begin
    isconst = T <: OrdinaryDiffEqConstantCache
    isconst || (cache = cache.const_cache)
    @unpack z,step = cache
    # This can be parallelized
    if !rewind
      if isconst
        for i in 1:step, j in step:-1:i
          z[j] = z[j] + z[j+1]
        end
      else
        for i in 1:step, j in step:-1:i
          @. z[j] = z[j] + z[j+1]
        end
      end # endif const cache
    else
      if isconst
        for i in 1:step, j in step:-1:i
          z[j] = z[j] - z[j+1]
        end
      else
        for i in 1:step, j in step:-1:i
          @. z[j] = z[j] - z[j+1]
        end
      end # endif const cache
    end # endif !rewind
  end # end @inbounds
end

# Apply corrections on the Nordsieck vector
function update_nordsieck_vector!(cache::T) where T
  isvode = ( T <: JVODECache || T <: JVODEConstantCache )
  ispreparevarorder = nordsieck_change_order(cache, 1)
  @inbounds begin
    isconst = T <: OrdinaryDiffEqConstantCache
    if isconst
      @unpack z,Δ,l,step = cache
      for i in 1:step+1
        z[i] = muladd.(l[i], Δ, z[i])
      end
      ispreparevarorder && ( z[end] = Δ )
    else
      @unpack z,Δ,l,step = cache.const_cache
      for i in 1:step+1
        @. z[i] = muladd(l[i], Δ, z[i])
      end
      ispreparevarorder && ( z[end] .= Δ )
    end # endif not const cache
  end # end @inbounds
end

function nlsolve_functional!(integrator, cache::T) where T
  @unpack f, dt, uprev, t, p = integrator
  isconstcache = T <: OrdinaryDiffEqConstantCache
  if isconstcache
    @unpack z, l, c_conv = cache
    ratetmp = integrator.f(z[1], p, dt+t)
  else
    @unpack ratetmp, const_cache = cache
    @unpack Δ, z, l, c_conv = const_cache
    cache = const_cache
    integrator.f(ratetmp, z[1], p, dt+t)
  end
  max_iter = 3
  div_rate = 2
  # Zero out the difference vector
  isconstcache ? ( cache.Δ = zero(cache.Δ) ) : ( Δ .= zero(eltype(Δ)) )
  # `k` is a counter for convergence test
  k = 0
  # `conv_rate` is used in convergence rate estimation
  conv_rate = 1.
  # initialize `δ_prev`
  δ_prev = 0
  # Start the functional iteration & store the difference into `Δ`
  for k in 1:max_iter
    if isconstcache
      ratetmp = inv(l[2])*muladd.(dt, ratetmp, -z[2])
      integrator.u = ratetmp + z[1]
      cache.Δ = ratetmp - cache.Δ
    else
      @. integrator.u = -z[2]
      @. ratetmp = inv(l[2])*muladd(dt, ratetmp, integrator.u)
      @. integrator.u = ratetmp + z[1]
      @. cache.Δ = ratetmp - cache.Δ
    end
    # @show norm(dt*ratetmp - ( z[2] + (integrator.u - z[1])*l[2] ))
    # @show norm(cache.Δ - (integrator.u - z[1]))
    # It only makes sense to calculate convergence rate in the second iteration
    δ = integrator.opts.internalnorm(cache.Δ)
    isconstcache ? ( cache.Δ = copy(ratetmp) ) : copy!(cache.Δ, ratetmp)
    if k >= 1
      conv_rate = max(1//10*conv_rate, δ/δ_prev)
    end
    test_rate = δ * min(one(conv_rate), conv_rate) / c_conv
    if test_rate <= one(test_rate)
      return true
    end
    # Divergence criteria
    if ( (k == max_iter) || (k >= 2 && δ > div_rate * δ_prev) )
      return false
    end
    δ_prev = δ
    isconstcache ? (ratetmp = integrator.f(integrator.u, p, dt+t)) :
                    integrator.f(ratetmp, integrator.u, p, dt+t)
  end
end

function nordsieck_rescale!(cache::T, rewind=false) where T
  isconstcache = T <: OrdinaryDiffEqConstantCache
  isconstcache || ( cache = cache.const_cache )
  @unpack z, tau, step = cache
  order = step
  eta = rewind ? tau[2]/tau[1] : tau[1]/tau[2]
  factor = eta
  for i in 2:order+1
    if isconstcache
      z[i] = z[i]*factor
    else
      scale!(z[i], factor)
    end
    factor *= eta
  end
  return nothing
end

function nordsieck_rewind!(cache)
  perform_predict!(cache, true)
  nordsieck_rescale!(cache, true)
end

function nordsieck_change_order(cache::T, n=0) where T
  isconstcache = T <: OrdinaryDiffEqConstantCache
  isconstcache || ( cache = cache.const_cache )
  isvode = ( T <: JVODECache || T <: JVODEConstantCache )
  isvode || return false
  cache.n_wait == 0+n
end

function nordsieck_decrement_wait!(cache::T) where T
  isvode = ( T <: JVODECache || T <: JVODEConstantCache )
  isvode || return nothing
  isconstcache = T <: OrdinaryDiffEqConstantCache
  isconstcache || ( cache = cache.const_cache )
  cache.n_wait = max(0, cache.n_wait-1)
  return nothing
end

function nordsieck_order_change(cache::T, dorder) where T
  isconstcache = T <: OrdinaryDiffEqConstantCache
  isconstcache || ( cache = cache.const_cache )
  @unpack step, tau = cache
  order = step
  # WIP: uncomment when finished
  #@inbound begin
  begin
    # Adams order increase
    if dorder == 1
      if isconstcache
        cache.z[order+2] = zero(cache.z[order+2])
      else
        cache.z[order+2] .= 0
      end
    else
      # Adams order decrease
      # One needs to rescale the Nordsieck vector on an order decrease
      cache.l .= 0
      cache.l[2] = 1
      dt = tau[1]
      hsum = zero(eltype(cache.tau))
      for j in 2:order-1
        hsum += cache.tau[j]
        # TODO: `hscale`?
        ξ = hsum / dt
        for i in j:-1:1
          cache.l[i+1] = cache.l[i+1] * ξ + cache.l[i]
        end # for i
      end # for j

      for j in 2:order-1
        cache.l[j+1] = order * cache.l[j] / j
      end
      for j in 3:order
        # cache.z[j] = -cache.l[j] * cache.z[order+1] + cache.z[j]
        if isconstcache
          cache.z[j] = muladd.(-cache.l[j], cache.z[order+1], cache.z[j])
        else
          @. cache.z[j] = muladd(-cache.l[j], cache.z[order+1], cache.z[j])
        end
      end # for j
    end # else
  end # @inbound
end

# `η` is `dtₙ₊₁/dtₙ`
function choose_η!(integrator, cache::T) where T
  isconstcache = T <: OrdinaryDiffEqConstantCache
  isconstcache || ( cache = cache.const_cache )
  isvarorder = nordsieck_change_order(cache)
  order = get_current_adaptive_order(integrator.alg, integrator.cache)
  L = order + 1
  ηq = stepsize_η!(integrator, cache, order)
  if isvarorder
    cache.n_wait = 2
    ηqm1 = stepsize_η₋₁!(integrator, cache, order)
    ηqp1 = stepsize_η₊₁!(integrator, cache, order)
    η = max(ηqm1, ηqp1, cache.η)
  else
    η = ηq
    cache.η = η
  end
  ( η <= integrator.opts.qsteady_max ) && ( cache.η = 1 ; return cache.η )
  if isvarorder
    if η == cache.η
      cache.nextorder = order
    elseif η == cache.η₋₁
      cache.η = cache.η₋₁
      cache.nextorder = order - 1
      cache.n_wait = L
      nordsieck_order_change(cache, -1)
    else
      cache.η = cache.η₊₁
      cache.nextorder = order + 1
      # TODO: BDF needs a different handler
      cache.n_wait = L
      nordsieck_order_change(cache, 1)
    end
  end
  ( integrator.iter == 1 || integrator.u_modified ) && return ( cache.η = min(1e5, cache.η) )
  cache.η = min(integrator.opts.qmax, max(integrator.opts.qmin, cache.η))
  return cache.η
end

function stepsize_η!(integrator, cache::T, order) where T
  L = order+1
  cache.η = inv( (BIAS2*integrator.EEst)^inv(L) + ADDON )
  return cache.η
end

function stepsize_η₊₁!(integrator, cache::T, order) where T
  isconstcache = T <: OrdinaryDiffEqConstantCache
  isconstcache || ( atmp = cache.atmp; cache = cache.const_cache )
  @unpack uprev, u = integrator
  @unpack z, c_LTE₊₁, tau, c_𝒟  = cache
  q = order
  cache.η₊₁ = 0
  qmax = length(z)-1
  L = q+1
  if q != qmax
    cache.prev_𝒟 == 0 && return cache.η₊₁
    cquot = (c_𝒟 / cache.prev_𝒟) * (tau[1]/tau[3])^L
    if isconstcache
      @show atmp = muladd.(-cquot, z[end], cache.Δ)
      @show atmp = calculate_residuals(atmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    else
      @. atmp = muladd(-cquot, z[end], cache.Δ)
      calculate_residuals!(atmp, const_cache.Δ, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    end
    @show dup = abs(integrator.opts.internalnorm(atmp) * c_LTE₊₁)
    cache.η₊₁ = inv( (BIAS3*dup)^inv(L+1) + ADDON )
  end
  return cache.η₊₁
end

function stepsize_η₋₁!(integrator, cache::T, order) where T
  isconstcache = T <: OrdinaryDiffEqConstantCache
  isconstcache || ( atmp = cache.atmp; cache = cache.const_cache )
  @unpack uprev, u = integrator
  @unpack z, c_LTE₋₁ = cache
  q = order
  cache.η₋₁ = 0
  if q > 1
    if isconstcache
      atmp = calculate_residuals(cache.z[q+1], uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    else
      calculate_residuals!(atmp, const_cache.Δ, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    end
    approx = integrator.opts.internalnorm(atmp) * c_LTE₋₁
    cache.η₋₁ = inv( (BIAS1*approx)^inv(q) + ADDON )
  end
  return cache.η₋₁
end
