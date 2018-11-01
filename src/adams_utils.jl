# Solving Ordinary Differential Equations I: Nonstiff Problems
# by Ernst Hairer, Gerhard Wanner, and Syvert P Norsett.
# III.5 Variable Step Size Multistep Methods: Formulae 5.9
function ϕ_and_ϕstar!(cache, du, k)
  @inbounds begin
    @unpack dts, ϕstar_nm1, ϕ_n, ϕstar_n, β = cache
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
  end # inbounds
  nothing
end

function ϕ_and_ϕstar!(cache::Union{VCABMConstantCache,VCABMCache}, du, k)
  @inbounds begin
    @unpack dts, ϕstar_nm1, ϕ_n, ϕstar_n, β = cache
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
    cache.ξ = ξ
    cache.ξ0 = ξ0
  end # inbounds
  nothing
end

function expand_ϕ_and_ϕstar!(cache, i)
  @unpack ξ, ξ0, β, dts, ϕstar_nm1, ϕ_n, ϕstar_n = cache
  ξ0 += dts[i]
  β[i] = β[i-1] * ξ/ξ0
  if typeof(cache) <: OrdinaryDiffEqMutableCache
    @. ϕ_n[i] = ϕ_n[i-1] - ϕstar_nm1[i-1]
    @. ϕstar_n[i] = β[i] * ϕ_n[i]
  else
    ϕ_n[i] = ϕ_n[i-1] - ϕstar_nm1[i-1]
    ϕstar_n[i] = β[i] * ϕ_n[i]
  end
  nothing
end

function ϕ_np1!(cache, du_np1, k)
  @inbounds begin
    @unpack ϕ_np1, ϕstar_n = cache
    for i = 1:k
      if i != 1
        if typeof(cache) <: OrdinaryDiffEqMutableCache
          @. ϕ_np1[i] = ϕ_np1[i-1] - ϕstar_n[i-1]
        else
          ϕ_np1[i] = ϕ_np1[i-1] - ϕstar_n[i-1]
        end
      else
        if typeof(cache) <: OrdinaryDiffEqMutableCache
          ϕ_np1[i] .= du_np1
        else
          ϕ_np1[i] = du_np1
        end
      end
    end
  end #inbounds
  nothing
end

# Solving Ordinary Differential Equations I: Nonstiff Problems
# by Ernst Hairer, Gerhard Wanner, and Syvert P Norsett.
# III.5 Variable Step Size Multistep Methods: Formulae 5.9 & 5.10
# ------------------------------------------------------------
# Note that `g` is scaled by `dt` in here
function g_coefs!(cache, k)
  @inbounds begin
    @unpack dts,c,g = cache
    ξ = dt = dts[1]
    for i = 1:k
      if i > 2
        ξ += dts[i-1]
      end
      for q = 1:k-(i-1)
        if i > 2
          c[i,q] =  muladd(-dt/ξ, c[i-1,q+1], c[i-1,q])
        elseif i == 1
          c[i,q] = inv(q)
        elseif i == 2
          c[i,q] = inv(q*(q+1))
        end
      end # q
      g[i] = c[i,1] * dt
    end # i
  end # inbounds
  nothing
end

# Coefficients for the implicit Adams methods
const γstar = @SVector [1,-1/2,-1/12,-1/24,-19/720,-3/160,-863/60480,-275/24192,-33953/3628800,-0.00789255,-0.00678585,-0.00592406,-0.00523669,-0.0046775,-0.00421495,-0.0038269]
