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
