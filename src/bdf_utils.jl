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
