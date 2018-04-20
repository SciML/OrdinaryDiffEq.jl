# This function computes the integral, from -1 to 0, of a polynomial
# `P(x)` from the coefficients of `P` with an offset `k`.
function ∫₋₁⁰dx(a, deg, k)
  @inbounds begin
    int = zero(eltype(a))
    sign = one(eltype(a))
    for i in 1:deg
      int += sign * a[i]/(i+k)
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
function calc_coeff!(cache)
  @inbounds begin
    @unpack m, l, tau = cache
    ZERO, ONE = zero(m[1]), one(m[1])
    dtsum = dt = tau[1]
    order = cache.step
    m[1] = ONE
    for i in 2:order+1
      m[i] = ZERO
    end
    for j in 1:order-1
      ξ_inv = dt / dtsum
      for i in j:-1:1
        m[i+1] += m[i] * ξ_inv
      end
      dtsum += tau[j+1]
    end

    M0 = ∫₋₁⁰dx(m, order, 0)
    M1 = ∫₋₁⁰dx(m, order, 1)
    M0_inv = inv(M0)
    l[1] = ONE
    for i in 1:order
      l[i+1] = M0_inv * m[i] / i
    end
    cache.tq = M1 * M0_inv * ξ_inv
  end
end

# Apply the Pascal linear operator
function perform_predict!(cache, undo)
  @inbounds begin
    @unpack z,step = cache
    # This can be parallelized
    if !undo
      for i in 1:step, j in step:-1:i
        @. z[j] = z[j] + z[j+1]
      end
    else
      for i in 1:step, j in step:-1:i
        @. z[j] = z[j] - z[j+1]
      end
    end
  end
end

# Apply corrections on the Nordsieck vector
function perform_correct!(cache)
  @inbounds begin
    @unpack z,Δ,l,step = cache
    for i in 1:step+1
      @. z[i] = muladd(l[i], Δ, z[i])
    end
  end
end

# TODO: Functional iteration solver
function nlsolve_functional!(integrator, cache)
  @unpack f, dt, u, t, p = integrator
  @unpack tmp, tab = cache
  @unpack Δ, z, l, tq = tab
  integrator.f(tmp, z[1], p, dt+t)
  # Zero out the difference vector
  Δ .= zero(eltype(Δ))
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
    @. tmp = inv(l[2])*muladd(dt, tmp, -z[2])
    @. u = tmp + z[1]
    @. Δ = tmp - Δ
    δ = integrator.opt.internalnorm(Δ)
    # It only make sense to calcluate convergence rate in the second iteration
    if k >= 1
      conv_rate = max(1//10*conv_rate, δ/δ_prev)
    end
    test_rate = δ * min(one(conv_rate), conv_rate) / pconv
    test_rate <= one(test_rate) && return nothing
    k += 1
    # TODO: Add divergence test & max_iter
    ######################################
    δ_prev = δ
    integrator.f(tmp, u, p, dt+t)
  end
end
