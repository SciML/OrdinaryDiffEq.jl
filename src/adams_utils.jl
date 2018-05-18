# Solving Ordinary Differential Equations I: Nonstiff Problems
# by Ernst Hairer, Gerhard Wanner, and Syvert P Norsett.
# III.5 Variable Step Size Multistep Methods: Formulae 5.9
function ϕ_and_ϕstar!(cache, du)
  @inbounds begin
    @unpack dts, ϕstar_nm1, ϕ_n, ϕstar_n,β,k = cache
    ξ = dt = dts[1]
    ξ0 = zero(dt)
    β[1] = one(dt)
    if typeof(cache) <: OrdinaryDiffEqMutableCache
      ϕ_n[1] .= du
      ϕstar_n[1] .= du
    else
      ϕ_n[1] = du
      ϕstar_n[1] = du
    end
    for i = 2:k
      ξ0 += dts[i]
      β[i] = β[i-1] * ξ/ξ0
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

# Solving Ordinary Differential Equations I: Nonstiff Problems
# by Ernst Hairer, Gerhard Wanner, and Syvert P Norsett.
# III.5 Variable Step Size Multistep Methods: Formulae 5.9 & 5.10
function g_coefs!(cache)
  @inbounds begin
    @unpack dts,c,g,k = cache
    ξ = dt = dts[1]
    for i = 1:k
      if i > 2
        ξ += dts[i-1]
      end
      for q = 1:k-(i-1)
        if i == 1
          c[i,q] = inv(q)
        elseif i == 2
          c[i,q] = inv(q*(q+1))
        else
          c[i,q] = c[i-1,q] - c[i-1,q+1] * dt/ξ
        end
      end # q
      g[i] = c[i,1]
    end # i
    return g
  end # inbounds
end
