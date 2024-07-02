@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{DPRKN6ConstantCache, DPRKN6Cache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
@dprkn6pre0
return ArrayPartition(
    duprev +
    dt * Θ *
    (bp1Θ * k1 + bp3Θ * k3 +
     bp4Θ * k4 + bp5Θ * k5 + bp6Θ * k6),
    uprev +
    dt * Θ *
    (duprev +
     dt * Θ * (b1Θ * k1 + b3Θ * k3 +
               b4Θ * k4 + b5Θ * k5 + b6Θ * k6)))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{DPRKN6ConstantCache, DPRKN6Cache}, idxs,
    T::Type{Val{0}}, differential_vars::Nothing)
@dprkn6pre0
return ArrayPartition(
    duprev[idxs] +
    dt * Θ *
    (bp1Θ * k1[idxs] + bp3Θ * k3[idxs] +
     bp4Θ * k4[idxs] + bp5Θ * k5[idxs] + bp6Θ * k6[idxs]),
    uprev[idxs] +
    dt * Θ *
    (duprev[idxs] +
     dt * Θ *
     (b1Θ * k1[idxs] +
      b3Θ * k3[idxs] +
      b4Θ * k4[idxs] + b5Θ * k5[idxs] + b6Θ * k6[idxs])))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{DPRKN6ConstantCache, DPRKN6Cache}, idxs::Number,
    T::Type{Val{0}}, differential_vars::Nothing)
@dprkn6pre0
halfsize = length(y₀) ÷ 2
if idxs <= halfsize
    duprev[idxs] +
    dt * Θ *
    (bp1Θ * k1[idxs] + bp3Θ * k3[idxs] +
     bp4Θ * k4[idxs] + bp5Θ * k5[idxs] + bp6Θ * k6[idxs])
else
    idxs = idxs - halfsize
    uprev[idxs] +
    dt * Θ *
    (duprev[idxs] +
     dt * Θ *
     (b1Θ * k1[idxs] +
      b3Θ * k3[idxs] +
      b4Θ * k4[idxs] + b5Θ * k5[idxs] + b6Θ * k6[idxs]))
end
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{DPRKN6ConstantCache, DPRKN6Cache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
@dprkn6pre0
@inbounds @.. broadcast=false out.x[2]=uprev +
                                       dt * Θ *
                                       (duprev +
                                        dt * Θ *
                                        (b1Θ * k1 +
                                         b3Θ * k3 +
                                         b4Θ * k4 + b5Θ * k5 + b6Θ * k6))
@inbounds @.. broadcast=false out.x[1]=duprev +
                                       dt * Θ *
                                       (bp1Θ * k1 + bp3Θ * k3 +
                                        bp4Θ * k4 + bp5Θ * k5 + bp6Θ * k6)
#for i in eachindex(out.x[1])
#  out.x[2][i]  = uprev[i] + dt*Θ*(duprev[i] + dt*Θ*(b1Θ*k1[i] +
#                                                    b3Θ*k3[i] +
#                                                    b4Θ*k4[i] + b5Θ*k5[i] + b6Θ*k6[i]))
#  out.x[1][i] =  duprev[i] + dt*Θ*(bp1Θ*k1[i] + bp3Θ*k3[i] +
#                                   bp4Θ*k4[i] + bp5Θ*k5[i] + bp6Θ*k6[i])
#end
out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{DPRKN6ConstantCache, DPRKN6Cache}, idxs,
    T::Type{Val{0}}, differential_vars::Nothing)
@dprkn6pre0
halfsize = length(y₀) ÷ 2
isfirsthalf = idxs .<= halfsize
secondhalf = idxs .> halfsize
firstidxs = idxs[isfirsthalf]
secondidxs_shifted = idxs[secondhalf]
secondidxs = secondidxs_shifted .- halfsize

@views @.. broadcast=false out[secondhalf]=uprev[secondidxs] +
                                           dt * Θ *
                                           (duprev[secondidxs] +
                                            dt * Θ *
                                            (b1Θ * k1[secondidxs] +
                                             b3Θ * k3[secondidxs] +
                                             b4Θ * k4[secondidxs] +
                                             b5Θ * k5[secondidxs] +
                                             b6Θ * k6[secondidxs]))
@views @.. broadcast=false out[isfirsthalf]=duprev[firstidxs] +
                                            dt * Θ *
                                            (bp1Θ * k1[firstidxs] +
                                             bp3Θ * k3[firstidxs] +
                                             bp4Θ * k4[firstidxs] +
                                             bp5Θ * k5[firstidxs] +
                                             bp6Θ * k6[firstidxs])
out
end