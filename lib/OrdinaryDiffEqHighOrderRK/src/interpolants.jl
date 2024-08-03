"""
"""
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP8ConstantCache, DP8Cache}, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing)
    Θ1 = 1 - Θ
    # return @.. broadcast=false y₀ + dt*Θ*(k[1] + Θ1*(k[2] + Θ*(k[3]+Θ1*(k[4] + Θ*(k[5] + Θ1*(k[6]+Θ*k[7]))))))
    return @inbounds y₀ +
                     dt * Θ *
                     (k[1] +
                      Θ1 * (k[2] +
                       Θ * (k[3] + Θ1 * (k[4] + Θ * (k[5] + Θ1 * (k[6] + Θ * k[7]))))))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP8ConstantCache, DP8Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing)
    Θ1 = 1 - Θ
    # return @.. broadcast=false y₀[idxs] + dt*Θ*(k[1][idxs] + Θ1*(k[2][idxs] + Θ*(k[3][idxs]+Θ1*(k[4][idxs] + Θ*(k[5][idxs] + Θ1*(k[6][idxs]+Θ*k[7][idxs]))))))
    return y₀[idxs] +
           dt * Θ *
           (k[1][idxs] +
            Θ1 * (k[2][idxs] +
             Θ * (k[3][idxs] +
              Θ1 * (k[4][idxs] + Θ * (k[5][idxs] + Θ1 * (k[6][idxs] + Θ * k[7][idxs]))))))
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP8ConstantCache, DP8Cache}, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing)
    Θ1 = 1 - Θ
    @inbounds @.. broadcast=false out=y₀ +
                                      dt * Θ *
                                      (k[1] +
                                       Θ1 * (k[2] +
                                        Θ * (k[3] +
                                         Θ1 *
                                         (k[4] + Θ * (k[5] + Θ1 * (k[6] + Θ * k[7]))))))
    #@inbounds for i in eachindex(out)
    #  out[i] = y₀[i] + dt*Θ*(k[1][i] + Θ1*(k[2][i] + Θ*(k[3][i]+Θ1*(k[4][i] + Θ*(k[5][i] + Θ1*(k[6][i]+Θ*k[7][i]))))))
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP8ConstantCache, DP8Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing)
    Θ1 = 1 - Θ
    @views @.. broadcast=false out=y₀[idxs] +
                                   dt * Θ *
                                   (k[1][idxs] +
                                    Θ1 * (k[2][idxs] +
                                     Θ * (k[3][idxs] +
                                      Θ1 * (k[4][idxs] +
                                       Θ *
                                       (k[5][idxs] + Θ1 * (k[6][idxs] + Θ * k[7][idxs]))))))
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = y₀[i] + dt*Θ*(k[1][i] + Θ1*(k[2][i] + Θ*(k[3][i]+Θ1*(k[4][i] + Θ*(k[5][i] + Θ1*(k[6][i]+Θ*k[7][i]))))))
    #end
    out
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP8ConstantCache, DP8Cache}, idxs::Nothing,
        T::Type{Val{1}}, differential_vars::Nothing)
    @inbounds b1diff = @.. broadcast=false k[1]+k[2]
    @inbounds b2diff = @.. broadcast=false -2*k[2]+2*k[3]+2*k[4]
    @inbounds b3diff = @.. broadcast=false -3 * k[3]-6 * k[4]+3*k[5]+3*k[6]
    @inbounds b4diff = @.. broadcast=false 4 * k[4] - 8 * k[5] - 12 * k[6]+4 * k[7]
    @inbounds b5diff = @.. broadcast=false 5 * k[5] + 15 * k[6]-15 * k[7]
    @inbounds b6diff = @.. broadcast=false -6 * k[6]+18 * k[7]
    @inbounds b7diff = @.. broadcast=false -7*k[7]
    # return @.. broadcast=false b1diff + Θ*(b2diff + Θ*(b3diff + Θ*(b4diff + Θ*(b5diff + Θ*(b6diff + Θ*b7diff)))))
    return b1diff +
           Θ *
           (b2diff + Θ * (b3diff + Θ * (b4diff + Θ * (b5diff + Θ * (b6diff + Θ * b7diff)))))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP8ConstantCache, DP8Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing)
    b1diff = @.. broadcast=false k[1][idxs]+k[2][idxs]
    b2diff = @.. broadcast=false -2*k[2][idxs]+2*k[3][idxs]+2*k[4][idxs]
    b3diff = @.. broadcast=false -3 * k[3][idxs]-6 * k[4][idxs]+3*k[5][idxs]+3*k[6][idxs]
    b4diff = @.. broadcast=false 4 * k[4][idxs] - 8 * k[5][idxs] -
                                 12 * k[6][idxs]+4 * k[7][idxs]
    b5diff = @.. broadcast=false 5 * k[5][idxs] + 15 * k[6][idxs]-15 * k[7][idxs]
    b6diff = @.. broadcast=false -6 * k[6][idxs]+18 * k[7][idxs]
    b7diff = @.. broadcast=false -7*k[7][idxs]
    # return @.. broadcast=false b1diff + Θ*(b2diff + Θ*(b3diff + Θ*(b4diff + Θ*(b5diff + Θ*(b6diff + Θ*b7diff)))))
    return b1diff +
           Θ *
           (b2diff + Θ * (b3diff + Θ * (b4diff + Θ * (b5diff + Θ * (b6diff + Θ * b7diff)))))
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP8ConstantCache, DP8Cache}, idxs::Nothing,
        T::Type{Val{1}}, differential_vars::Nothing)
    # b1diff = k[1] + k[2]
    # b2diff = -2*k[2] + 2*k[3] + 2*k[4]
    # b3diff = -3*k[3] - 6*k[4] + 3*k[5] + 3*k[6]
    # b4diff = 4*k[4] - 8*k[5] - 12*k[6] + 4*k[7]
    # b5diff = 5*k[5] + 15*k[6] - 15*k[7]
    # b6diff = -6*k[6] + 18*k[7]
    # @.. broadcast=false out = b1diff + Θ*(b2diff + Θ*(b3diff + Θ*(b4diff +
    #                                                Θ*(b5diff + Θ*(b6diff - 7*k[7]*Θ)))))
    @views @.. broadcast=false out=k[1] + k[2] +
                                   Θ * (-2 * k[2] + 2 * k[3] + 2 * k[4] +
                                    Θ * (-3 * k[3] - 6 * k[4] + 3 * k[5] + 3 * k[6] +
                                     Θ * (4 * k[4] - 8 * k[5] - 12 * k[6] + 4 * k[7] +
                                      Θ * (5 * k[5] + 15 * k[6] - 15 * k[7] +
                                       Θ * (-6 * k[6] + 18 * k[7] - 7 * k[7] * Θ)))))
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP8ConstantCache, DP8Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing)
    # b1diff = k[1][idxs] + k[2][idxs]
    # b2diff = -2*k[2][idxs] + 2*k[3][idxs] + 2*k[4][idxs]
    # b3diff = -3*k[3][idxs] - 6*k[4][idxs] + 3*k[5][idxs] + 3*k[6][idxs]
    # b4diff = 4*k[4][idxs] - 8*k[5][idxs] - 12*k[6][idxs] + 4*k[7][idxs]
    # b5diff = 5*k[5][idxs] + 15*k[6][idxs] - 15*k[7][idxs]
    # b6diff = -6*k[6][idxs] + 18*k[7][idxs]
    #@views @.. broadcast=false out = b1diff + Θ*(b2diff + Θ*(b3diff + Θ*(b4diff +
    #                                               Θ*(b5diff + Θ*(b6diff - 7*k[7][idxs]*Θ)))))
    @views @.. broadcast=false out=k[1][idxs] + k[2][idxs] +
                                   Θ * (-2 * k[2][idxs] + 2 * k[3][idxs] + 2 * k[4][idxs] +
                                    Θ *
                                    (-3 * k[3][idxs] - 6 * k[4][idxs] + 3 * k[5][idxs] +
                                     3 * k[6][idxs] +
                                     Θ *
                                     (4 * k[4][idxs] - 8 * k[5][idxs] - 12 * k[6][idxs] +
                                      4 * k[7][idxs] +
                                      Θ *
                                      (5 * k[5][idxs] + 15 * k[6][idxs] - 15 * k[7][idxs] +
                                       Θ * (-6 * k[6][idxs] + 18 * k[7][idxs] -
                                        7 * k[7][idxs] * Θ)))))
    out
end