@inline function DiffEqBase.get_du!(out, integrator::ODEIntegrator)
    integrator.cache isa FunctionMapCache ||
        integrator.cache isa FunctionMapConstantCache &&
            error("Derivatives are not defined for this stepper.")
    if integrator.cache isa FunctionMapCache
        out .= integrator.cache.tmp
    else
        return if isdefined(integrator, :fsallast) &&
                  !(integrator.alg isa
                    Union{Rosenbrock23, Rosenbrock32, Rodas23W,
            Rodas3P, Rodas4, Rodas4P, Rodas4P2, Rodas5,
            Rodas5P, Rodas5Pe, Rodas5Pr})
            # Special stiff interpolations do not store the right value in fsallast
            out .= integrator.fsallast
        else
            integrator(out, integrator.t, Val{1})
        end
    end
end

@inline function DiffEqBase.get_tmp_cache(integrator,
        alg::OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
        cache::OrdinaryDiffEqMutableCache)
    (cache.tmp, cache.linsolve_tmp)
end

function resize_non_user_cache!(integrator::ODEIntegrator,
        cache::RosenbrockMutableCache, i)
    cache.J = similar(cache.J, i, i)
    cache.W = similar(cache.W, i, i)
    resize_jac_config!(cache.jac_config, i)
    resize_grad_config!(cache.grad_config, i)
    nothing
end