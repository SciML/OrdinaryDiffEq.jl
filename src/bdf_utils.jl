# bdf_utils

function U!(k, U, inplace)
  for j = 1:k
    for r = 1:k
      if inplace
        @. U[j,r] = inv(factorial(j)) * prod([m-r for m in 0:(j-1)])
      else
        U[j,r] = inv(factorial(j)) * prod([m-r for m in 0:(j-1)]) 
      end
    end
  end  
end

function R!(k, ρ, cache)
  @unpack R = cache
  for j = 1:k
    for r = 1:k
      if typeof(cache) <: OrdinaryDiffEqMutableCache
        @. R[j,r] = inv(factorial(j)) * prod([m-r*ρ for m in 0:(j-1)])
      else
        R[j,r] = inv(factorial(j)) * prod([m-r*ρ for m in 0:(j-1)])
      end
    end
  end
end
