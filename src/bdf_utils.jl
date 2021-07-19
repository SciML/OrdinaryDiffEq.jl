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

error_constant(integrator, order) = error_constant(integrator, integrator.alg, order)
function error_constant(integrator, alg::QNDF, k)
  @unpack γₖ = integrator.cache
  κ = alg.kappa[k]
  κ*γₖ[k] + inv(k+1)
end

function compute_weights!(ts,k,weights)
  for i = 1:k+1
      weights[i] = one(eltype(weights))
      for j = 1:i-1
          weights[i] *= ts[i] - ts[j]
      end
      for j = i+1:k+1
          weights[i] *= ts[i] - ts[j]
      end
      weights[i] = 1/weights[i]
  end
end

function calc_Lagrange_interp(k,weights,t,ts,u_history,u::Number)
  #@show t,ts,u_history
  if t in ts
    i = searchsortedfirst(ts,t,rev=true)
    return u_history[i]
  else
    for i in 1:k+1
      u += weights[i]/(t-ts[i])*u_history[i]
    end
    for i in 1:k+1
      u *= t-ts[i]
    end
  end
  #@show weights
  u
end

function calc_Lagrange_interp(k,weights,t,ts,u_history,u)
  if t in ts
    i = searchsortedfirst(ts,t,rev=true)
    return u_history[:,i]
  else
    for i in 1:k+1
      @.. u += weights[i]/(t-ts[i])*u_history[:,i]
    end
    for i in 1:k+1
      @.. u *= t-ts[i]
    end
  end
  u
end

function calc_Lagrange_interp!(k,weights,t,ts,u_history,u)
  if t in ts
    i = searchsortedfirst(ts,t,rev=true)
    @.. u = u_history[:,i]
  else
    for i in 1:k+1
      @.. u += weights[i]/(t-ts[i])*u_history[:,i]
    end
    for i in 1:k+1
      @.. u *= t-ts[i]
    end
  end
end

#This code refers to https://epubs.siam.org/doi/abs/10.1137/S0036144596322507

function calc_finite_difference_weights(ts,t,order,::Val{N}) where {N}
  max_order = N
  c = zero(MMatrix{max_order+1,max_order+1,eltype(ts)})
  c1 = one(t)
  c4 = ts[1] - t
  c[1,1] = one(t)
  for i in 2:order+1
    c2 = one(t)
    c5 = c4
    c4 = ts[i] - t
    @inbounds for j = 1:i-1
      c3 = ts[i] - ts[j]
      c2 *= c3
      if j == i-1
        for k in i:-1:2
          c[i,k] = c1*((k-1)*c[i-1,k-1]-c5*c[i-1,k])/c2
        end
        c[i,1] = zero(t)
      end
      for k in i:-1:2
        c[j,k] = (c4*c[j,k] - (k-1)*c[j,k-1])/c3
      end
      c[j,1] = zero(t)
    end
    c1 = c2
  end
  return SArray(c)
end
