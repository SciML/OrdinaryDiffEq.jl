# bdf_utils
function U!(k, U)
  for j = 1:k
    for r = 1:k
      if j == 1
        U[j,r] = -r
      else
        U[j,r] = U[j-1,r] * ((j-1) - r)/j
      end
    end
  end
  nothing
end

function R!(k, ρ, cache)
  @unpack R = cache
  for j = 1:k
    for r = 1:k
      if j == 1
        R[j,r] = -r * ρ
      else
        R[j,r] = R[j-1,r] * ((j-1) - r * ρ)/j
      end
    end
  end
  nothing
end

# This functions takes help of D2 array to create backward differences array D
# Ith row of D2 keeps Ith order backward differences (∇ⁱyₙ)
function backward_diff!(cache, D, D2, k, flag=true)
  if flag
    if typeof(cache) <: OrdinaryDiffEqMutableCache
      D[1] .= D2[1,1]
    else
      D[1] = D2[1,1]
    end
  end
  for i = 2:k
    for j = 1:(k-i+1)
      if typeof(cache) <: OrdinaryDiffEqMutableCache
        @. D2[i,j] = D2[i-1,j] - D2[i-1,j+1]
      else
        D2[i,j] = D2[i-1,j] - D2[i-1,j+1]
      end
    end
    if flag
      if typeof(cache) <: OrdinaryDiffEqMutableCache
        D[i] .= D2[i,1]
      else
        D[i] = D2[i,1]
      end
    end
  end
  nothing
end

# this function updates backward difference array D when stepsize gets change
# Formula -> D = D * (R * U)
# and it is taken from the paper -
# Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s paper
# E. Alberdi Celaya, J. J. Anza Aguirrezabala, and P. Chatzipantelidis
function reinterpolate_history!(cache, D, R, k)
  @unpack tmp = cache
  if typeof(cache) <: OrdinaryDiffEqMutableCache
    fill!(tmp,zero(eltype(D[1])))
  else
    tmp = zero(eltype(D))
  end
  for j = 1:k
    for k = 1:k
      if typeof(cache) <: OrdinaryDiffEqMutableCache
        @. tmp += D[k] * R[k,j]
      else
        tmp += D[k] * R[k,j]
      end
    end
    if typeof(cache) <: OrdinaryDiffEqMutableCache
      D[j] .= tmp
      fill!(tmp,zero(eltype(D[1])))
    else
      D[j] = tmp
      tmp = zero(eltype(D))
    end
  end
  nothing
end

const γₖ = @SVector[sum(1//j for j in 1:k) for k in 1:6]

# this stepsize and order controller is taken from
# Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s paper
# E. Alberdi Celaya, J. J. Anza Aguirrezabala, and P. Chatzipantelidis
function stepsize_and_order!(cache, est, estₖ₋₁, estₖ₊₁, h, k)
  zₛ = 1.2
  zᵤ = 0.1
  Fᵤ = 10

  expo = 1/(k+1)
  z = zₛ * ((est)^expo)
  F = inv(z)

  hₖ₋₁ = 0.0
  hₖ₊₁ = 0.0

  if z <= zₛ
    # step is successful
    # precalculations
    if z <= zᵤ
      hₖ = Fᵤ * h
    elseif zᵤ < z <= zₛ
      hₖ = F * h
    end

    if k > 1
      expo = 1/k
      zₖ₋₁ = 1.3 * ((estₖ₋₁)^expo)
      Fₖ₋₁ = inv(zₖ₋₁)
      if zₖ₋₁ <= 0.1
        hₖ₋₁ = 10 * h
      elseif 0.1 < zₖ₋₁ <= 1.3
        hₖ₋₁ = Fₖ₋₁ * h
      end
    end

    expo = 1/(k+2)
    zₖ₊₁ = 1.4 * ((estₖ₊₁)^expo)
    Fₖ₊₁ = inv(zₖ₊₁)

    if zₖ₊₁<= 0.1
      hₖ₊₁ = 10 * h
    elseif 0.1 < zₖ₊₁ <= 1.4
      hₖ₊₁ = Fₖ₊₁ * h
    end
    # adp order and step conditions
    if hₖ₋₁ > hₖ
      hₙ = hₖ₋₁
      kₙ = max(k-1,1)
    else
      hₙ = hₖ
      kₙ = k
    end
    if hₖ₊₁ > hₙ
      hₙ = hₖ₊₁
      kₙ = min(k+1,5)
    end
    if hₙ < h
      hₙ = h
      kₙ = k
    end
    cache.h = hₙ
    cache.order = kₙ
    return true
  else
    # step is not successful
    if cache.c >= 1  # postfail
      cache.h = h/2
      cache.order = k
      return
    end
    if 1.2 < z <= 10
      hₖ = F * h
    elseif z > 10
      hₖ = 0.1 * h
    end
    hₙ = hₖ
    kₙ = k
    if k > 1
      expo = 1/k
      zₖ₋₁ = 1.3 * ((estₖ₋₁)^expo)
      Fₖ₋₁ = inv(zₖ₋₁)
      if 1.3 < zₖ₋₁ <= 10
        hₖ₋₁ = Fₖ₋₁ * h
      elseif zₖ₋₁ > 10
        hₖ₋₁ = 0.1 * h
      end

      if hₖ₋₁ > hₖ
        hₙ = min(h,hₖ₋₁)
        kₙ = max(k-1,1)
      end
    end
    cache.h = hₙ
    cache.order = kₙ
    return false
  end
end
