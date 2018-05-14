function ϕ_and_ϕstar!(cache, dy, next_point, last_idx)
  @unpack grid_points, ϕstar_nm1, ϕ_n, ϕstar_n,β,k = cache
  for i = 1:k
    if i == 1
      β[1] = 1
      if typeof(cache) <: OrdinaryDiffEqMutableCache
        ϕ_n[1] .= dy
        ϕstar_n[1] .= dy
      else
        ϕ_n[1] = dy
        ϕstar_n[1] = dy
      end
    else
      β[i] = β[i-1] * (next_point - grid_points[last_idx-i+2])/(grid_points[last_idx] - grid_points[last_idx-i+1])
      if typeof(cache) <: OrdinaryDiffEqMutableCache
        @. ϕ_n[i] = ϕ_n[i-1] - ϕstar_nm1[i-1]
        @. ϕstar_n[i] = β[i] * ϕ_n[i]
      else
        ϕ_n[i] = ϕ_n[i-1] - ϕstar_nm1[i-1]
        ϕstar_n[i] = β[i] * ϕ_n[i]
      end
    end
  end
  return ϕ_n, ϕstar_n
end

function g_coefs!(cache, dt, next_point, last_idx)
  @unpack grid_points,c,g,k = cache
  for i = 1:k
    for q = 1:(k-(i-1))
      if i == 1
        c[i,q] = 1/q
      elseif i == 2
        c[i,q] = 1/q/(q+1)
      else
        c[i,q] = c[i-1,q] - c[i-1,q+1] * dt/(next_point - grid_points[last_idx-i+1+1])
      end
    end
    g[i] = c[i,1]
  end
  return g
end
