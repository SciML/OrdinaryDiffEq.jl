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
end

# This functions takes help of D2 array to create backward differences array D
# Ith row of D2 keeps Ith order backward differences (∇ⁱyₙ)
function backward_diff!(D, D2, k, flag=true)
  if flag
    D[1] = D2[1,1]
  end
  for i = 2:k
    for j = 1:(k-i+1)
      D2[i,j] = D2[i-1,j] - D2[i-1,j+1]
    end
    if flag
      D[i] = D2[i,1]
    end
  end
end

global const γₖ = [1//1, 3//2, 11//6, 25//12, 137//60, 49//20]

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
    cache.k = kₙ
    return true
  else
    # step is not successful
    if cache.c >= 1  # postfail
      cache.h = h/2
      cache.k = k
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
    cache.k = kₙ
    return false
  end
end
