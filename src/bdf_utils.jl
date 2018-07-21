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

function backward_diff(udiff, D, D2, k)
  for i = 1:k
    D2[1,i] = udiff[i]
  end
  D[1] = D2[1,1]
  for i = 2:k
    for j = 1:(k-i+1)
      D2[i,j] = D2[i-1,j] - D2[i-1,j+1]
    end
    D[i] = D2[i,1]
  end
end

global const γₖ = [1//1, 3//2, 11//6, 25//12, 137//60, 49//20]
