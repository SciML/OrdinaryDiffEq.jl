@inline function DiffEqBase.get_du(integrator::ODEIntegrator)
    integrator.cache isa FunctionMapCache ||
        integrator.cache isa FunctionMapConstantCache &&
            error("Derivatives are not defined for this stepper.")
    return if isdefined(integrator, :fsallast)
        integrator.fsallast
    else
        integrator(integrator.t, Val{1})
    end
end

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