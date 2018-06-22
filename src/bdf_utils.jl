# bdf_utils

function U!(k, cache)
  @unpack U = cache
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

function prev_u(uprev, k, t, dt, cache)
  @unpack D = cache
  uprev2 = uprev + (D[1,1]) * -1
  uprev2
end

function D2!(u, uprev, k, cache)
  @unpack D, D2 = cache
  D2[1,1] = D[1,1]
  D2[1,2] = u - uprev
  D2[2,1] = D2[1,2] - D2[1,1]  
end
