# bdf_utils

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
  @views @.. D[:,k+2] = dd - D[:,k+1]
  @views @.. D[:,k+1] = dd
  for i in k:-1:1
    @views @.. D[:,i] = D[:,i] + D[:,i+1]
  end
  return nothing
end