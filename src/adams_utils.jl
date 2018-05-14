function ϕ_and_ϕstar!(cache, du)
  begin
    @unpack dts, ϕstar_nm1, ϕ_n, ϕstar_n,β,k = cache
    ξ = dt = dts[1]
    β[1] = 1
    if typeof(cache) <: OrdinaryDiffEqMutableCache
      ϕ_n[1] .= du
      ϕstar_n[1] .= du
    else
      ϕ_n[1] = du
      ϕstar_n[1] = du
    end
    for i = 2:k
      β[i] = β[i-1] * (dt+ξ)/ξ
      ξ += dts[i]
      if typeof(cache) <: OrdinaryDiffEqMutableCache
        @. ϕ_n[i] = ϕ_n[i-1] - ϕstar_nm1[i-1]
        @. ϕstar_n[i] = β[i] * ϕ_n[i]
      else
        ϕ_n[i] = ϕ_n[i-1] - ϕstar_nm1[i-1]
        ϕstar_n[i] = β[i] * ϕ_n[i]
      end
    end
    return ϕ_n, ϕstar_n
  end # inbounds
end

function g_coefs!(cache)
  begin
    @unpack dts,c,g,k = cache
    dt = dts[1]
    ξ = dt+dts[2]
    for i = 1:k
      for q = 1:k-(i-1)
        if i == 1
          c[i,q] = 1/q
        elseif i == 2
          c[i,q] = 1/q/(q+1)
        else
          c[i,q] = c[i-1,q] - c[i-1,q+1] * dt/ξ
          ξ += dts[i]
        end
      end
      g[i] = c[i,1]
    end
    return g
  end # inbounds
end
