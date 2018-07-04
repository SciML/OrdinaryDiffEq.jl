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

function prev_u!(uprev, k, t, dt, cache)
  @unpack D, uprev2 = cache
  if typeof(cache) <: OrdinaryDiffEqMutableCache
    @. uprev2 = uprev + (D[1,1]) * -1
  else
    uprev2 = uprev + (D[1,1]) * -1
    return uprev2
  end
end

function D2!(u, uprev, k, cache)
  @unpack D, D2 = cache
  if typeof(cache) <: OrdinaryDiffEqMutableCache
    @. D2[1,1] = (u - uprev) - D[1,1] 
  else
    D2[1,1] = (u - uprev) - D[1,1] 
  end
end
