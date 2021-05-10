# bdf_utils

function R!(k, ρ, cache)
  @unpack R = cache
  @inbounds for r = 1:k
    R[1,r] = -r * ρ
    for j = 2:k
      R[j,r] = R[j-1,r] * ((j-1) - r * ρ)/j
    end
  end
  nothing
end

# This functions takes help of D2 array to create backward differences array D
# Ith row of D2 keeps Ith order backward differences (∇ⁱyₙ)
function backward_diff!(cache::OrdinaryDiffEqMutableCache, D, D2, k, flag=true)
  flag && copyto!(D[1], D2[1,1])
  for i = 2:k
    for j = 1:(k-i+1)
      @.. D2[i,j] = D2[i-1,j] - D2[i-1,j+1]
    end
    flag && copyto!(D[i], D2[i,1])
  end
end

function backward_diff!(cache::OrdinaryDiffEqConstantCache, D, D2, k, flag=true)
  flag && (D[1] = D2[1,1])
  for i = 2:k
    for j = 1:(k-i+1)
      D2[i,j] = D2[i-1,j] - D2[i-1,j+1]
    end
    flag && (D[i] = D2[i,1])
  end
end

# this function updates backward difference array D when stepsize gets change
# Formula -> D = D * (R * U)
# and it is taken from the paper -
# Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s paper
# E. Alberdi Celaya, J. J. Anza Aguirrezabala, and P. Chatzipantelidis
function reinterpolate_history!(cache::OrdinaryDiffEqMutableCache, D, R, k)
  @unpack tmp = cache.nlsolver
  fill!(tmp,zero(eltype(D[1])))
  for j = 1:k
    for k = 1:k
      @. tmp += D[k] * R[k,j]
    end
    D[j] .= tmp
    fill!(tmp, zero(eltype(tmp)))
  end
end

function reinterpolate_history!(cache::OrdinaryDiffEqConstantCache, D, R, k)
  tmp = zero(D[1])
  for j = 1:k
    for k = 1:k
      tmp += D[k] * R[k,j]
    end
    D[j] = tmp
  end
end
