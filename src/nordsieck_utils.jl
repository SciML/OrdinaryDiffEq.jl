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
