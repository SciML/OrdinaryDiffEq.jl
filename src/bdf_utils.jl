# bdf_utils
@inline function U!(k, U)
  @inbounds for r = 1:k
    U[1,r] = -r
    for j = 2:k
      U[j,r] = U[j-1,r] * ((j-1) - r)/j
    end
  end
  nothing
end
  
function R!(k, ρ, cache)
  @unpack R = cache
  @inbounds for r = 1:k
    R[1,r] = -r * ρ
    for j = 2:k
      R[j,r] = R[j-1,r] * ((j-1) - r * ρ)/j
    end
  end
  nothing
end

# This functions takes help of D2 array to create backward differences array D
# Ith row of D2 keeps Ith order backward differences (∇ⁱyₙ)

function backward_diff!(cache::OrdinaryDiffEqMutableCache, D, D2, k, flag=true)
  flag && copyto!(D[1], D2[1,1])
  for i = 2:k
    for j = 1:(k-i+1)
      @.. D2[i,j] = D2[i-1,j] - D2[i-1,j+1]
    end
    flag && copyto!(D[i], D2[i,1])
  end
end

function backward_diff!(cache::OrdinaryDiffEqConstantCache, D, D2, k, flag=true)
  flag && (D[1] = D2[1,1])
  for i = 2:k
    for j = 1:(k-i+1)
      D2[i,j] = D2[i-1,j] - D2[i-1,j+1]
    end
    flag && (D[i] = D2[i,1])
  end
end

# this function updates backward difference array D when stepsize gets change
# Formula -> D = D * (R * U)
# and it is taken from the paper -
# Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s paper
# E. Alberdi Celaya, J. J. Anza Aguirrezabala, and P. Chatzipantelidis
function reinterpolate_history!(cache::OrdinaryDiffEqMutableCache, D, R, k)
  @unpack tmp = cache.nlsolver
  fill!(tmp,zero(eltype(D[1])))
  for j = 1:k
    for k = 1:k
      @. tmp += D[k] * R[k,j]
    end
    D[j] .= tmp
    fill!(tmp, zero(eltype(tmp)))
  end
end

function reinterpolate_history!(cache::OrdinaryDiffEqConstantCache, D, R, k)
  tmp = zero(D[1])
  for j = 1:k
    for k = 1:k
      tmp += D[k] * R[k,j]
    end
    D[j] = tmp
  end
end


function calc_R(ρ, k, ::Val{N}) where {N}
  R = zero(MMatrix{N,N,typeof(ρ)})
  @inbounds for r = 1:k
    R[1,r] = -r * ρ
    for j = 2:k
      R[j,r] = R[j-1,r] * ((j-1) - r * ρ)/j
    end
  end
  SArray(R)
end

function update_D!(D, dd, k)
  dd = _vec(dd)
  @views @.. D[:,k+2] = dd - D[:,k+1]
  @views @.. D[:,k+1] = dd
  for i in k:-1:1
    @views @.. D[:,i] = D[:,i] + D[:,i+1]
  end
  return nothing
end

const γₖ = @SVector[sum(1//j for j in 1:k) for k in 1:6]

###FBDF

function compute_weights(ts, y_history, ::Val{N}) where {N}
  ω = zero(MVector{N,eltype(y_history)})
  ω[1] = 1
  @inbounds for j in 2:N
    for k in 1:j-2
      ω[k] = (ts[k] - ts[j])ω[k]
    end
    ω[j] = 1
    for k in 1:j
      ω[j] *= (ts[j] - ts[k])
    end
  end
  @inbounds for j = 1:N
    ω[j] = inv(ω[j])
  end
  SArray(ω)
end

function interpolation_sols(weights, u::Vector, u_history, t, ts)
  for i in 1:k+1
    u += weights[i]/(t-ts[i])*u_history[i,:]
  end
  for i in 1:k+1
    u *= t+dt - ts[i]
  end
  u
end

function interpolation_sols!(weights, u::Matrix, u_history, t::Vector, ts)
  for j in 1:k
    for i in 1:k+1
      @.. u[j,:] += weights[i]/(t[j]-ts[i])*u_history[i,:]
    end
    for i in 1:k+1
      @.. u[j,:] *= t[j] - ts[i]
    end
  end
end
