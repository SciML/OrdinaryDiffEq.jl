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
  @inbounds begin
    isconst = T <: OrdinaryDiffEqConstantCache
    isconst || (cache = cache.const_cache)
    @unpack m, l, tau = cache
    ZERO, ONE = zero(m[1]), one(m[1])
    dtsum = dt = tau[1]
    order = cache.step
    m[1] = ONE
    for i in 2:order+1
      m[i] = ZERO
    end
    # initialize ξ_inv
    ξ_inv = dt / dtsum
    # compute coefficients from the Newton polynomial
    # check the `JuliaDiffEq/DiffEqDevMaterials` repository for more details
    for j in 1:order-1
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
    l[1] = ONE
    for i in 1:order
      l[i+1] = M0_inv * m[i] / i
    end
    cache.tq = order * l[order] * M1 * M0_inv * ξ_inv
  end
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
function perform_correct!(cache::T) where T
  @inbounds begin
    isconst = T <: OrdinaryDiffEqConstantCache
    if isconst
      @unpack z,Δ,l,step = cache
      for i in 1:step+1
        z[i] = muladd.(l[i], Δ, z[i])
      end
    else
      @unpack z,Δ,l,step = cache.const_cache
      for i in 1:step+1
        @. z[i] = muladd(l[i], Δ, z[i])
      end
    end # endif not const cache
  end # end @inbounds
end

function nlsolve_functional!(integrator, cache::T) where T
  @unpack f, dt, uprev, t, p = integrator
  isconstcache = T <: OrdinaryDiffEqConstantCache
  if isconstcache
    @unpack z, l, tq = cache
    ratetmp = integrator.f(z[1], p, dt+t)
  else
    @unpack ratetmp, const_cache = cache
    @unpack Δ, z, l, tq = const_cache
    cache = const_cache
    integrator.f(ratetmp, z[1], p, dt+t)
  end
  max_iter = 3
  div_rate = 2
  # Zero out the difference vector
  isconstcache ? ( cache.Δ = zero(cache.Δ) ) : ( Δ .= zero(eltype(Δ)) )
  # `pconv` is used in the convergence test
  pconv = (1//10) / tq
  # `k` is a counter for convergence test
  k = 0
  # `conv_rate` is used in convergence rate estimation
  conv_rate = 1.
  # initialize `δ_prev`
  δ_prev = 0
  # Start the functional iteration & store the difference into `Δ`
  while true
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
    k == 0 || isconstcache ? ( cache.Δ = copy(ratetmp) ) : copy!(cache.Δ, ratetmp)
    # It only makes sense to calculate convergence rate in the second iteration
    δ = integrator.opts.internalnorm(cache.Δ)
    if k >= 1
      conv_rate = max(1//10*conv_rate, δ/δ_prev)
    end
    test_rate = δ * min(one(conv_rate), conv_rate) / pconv
    test_rate <= one(test_rate) && return true
    k += 1
    # Divergence criteria
    (k == max_iter) || (k >= 2 && δ > div_rate * δ_prev) && return false
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
  eta = rewind ? tau[1]/tau[2] : tau[1]/tau[2]
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
