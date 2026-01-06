@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{FunctionMapConstantCache, FunctionMapCache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    y₀
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{FunctionMapConstantCache, FunctionMapCache},
        idxs, T::Type{Val{0}}, differential_vars::Nothing
    )
    y₀[idxs]
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{FunctionMapConstantCache, FunctionMapCache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    recursivecopy!(out, y₀)
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{FunctionMapConstantCache, FunctionMapCache},
        idxs, T::Type{Val{0}}, differential_vars::Nothing
    )
    @views out[idxs] .= y₀[idxs]
end
