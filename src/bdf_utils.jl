# bdf_utils

function U!(k, U)
  for j = 1:k
    for r = 1:k
      U[j,r] = inv(factorial(j)) * prod([m-r for m in 0:(j-1)]) 
    end
  end  
end

function R!(k, ρ, cache)
  @unpack R = cache
  for j = 1:k
    for r = 1:k
      R[j,r] = inv(factorial(j)) * prod([m-r*ρ for m in 0:(j-1)])
    end
  end
end
