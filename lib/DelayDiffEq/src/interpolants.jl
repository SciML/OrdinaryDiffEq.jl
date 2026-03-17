"""
    constant_extrapolant(t, integrator::DEIntegrator, idxs, deriv)

Calculate constant extrapolation of derivative `deriv` at time `t` and indices `idxs`
for `integrator`.
"""
function constant_extrapolant(t, integrator::DEIntegrator, idxs, deriv)
    return [constant_extrapolant(τ, integrator, idxs, deriv) for τ in t]
end

function constant_extrapolant(t::Number, integrator::DEIntegrator, idxs, T::Type{Val{0}})
    return idxs === nothing ? integrator.u : integrator.u[idxs]
end

function constant_extrapolant(t::Number, integrator::DEIntegrator, idxs, T::Type{Val{1}})
    return idxs === nothing ? zero(integrator.u) ./ oneunit(t) :
        zero(integrator.u[idxs]) ./ oneunit(t)
end

"""
    constant_extrapolant!(val, t, integrator::DEIntegrator, idxs, deriv)

Calculate constant extrapolation of derivative `deriv` at time `t` and indices `idxs`
for `integrator`, and save result in `val` if `val` is not `nothing`.
"""
function constant_extrapolant!(val, t, integrator::DEIntegrator, idxs, deriv)
    return [constant_extrapolant!(val, τ, integrator, idxs, deriv) for τ in t]
end

function constant_extrapolant!(
        val, t::Number, integrator::DEIntegrator, idxs,
        T::Type{Val{0}}
    )
    return if val === nothing
        if idxs === nothing
            return integrator.u
        else
            return integrator.u[idxs]
        end
    elseif idxs === nothing
        @. val = integrator.u
    else
        @views @. val = integrator.u[idxs]
    end
end

function constant_extrapolant!(
        val, t::Number, integrator::DEIntegrator, idxs,
        T::Type{Val{1}}
    )
    return if val === nothing
        if idxs === nothing
            return zero(integrator.u) ./ oneunit(t)
        else
            return zero(integrator.u[idxs]) ./ oneunit(t)
        end
    else
        fill!(val, zero(eltype(integrator.u)) / oneunit(t))
    end
end
