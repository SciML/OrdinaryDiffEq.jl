function DiffEqBase.initialize_dae!(
        integrator::Union{SDEIntegrator, AbstractSDDEIntegrator},
        initializealg = integrator.initializealg
    )
    return OrdinaryDiffEqCore._initialize_dae!(
        integrator, integrator.sol.prob, initializealg,
        Val(DiffEqBase.isinplace(integrator.sol.prob))
    )
end

function OrdinaryDiffEqCore._initialize_dae!(
        integrator::Union{SDEIntegrator, AbstractSDDEIntegrator},
        prob::Union{SciMLBase.AbstractRODEProblem, SciMLBase.AbstractSDDEProblem},
        ::OrdinaryDiffEqCore.DefaultInit, isinplace
    )
    return if SciMLBase.has_initializeprob(prob.f)
        OrdinaryDiffEqCore._initialize_dae!(integrator, prob, SciMLBase.OverrideInit(), isinplace)
    elseif SciMLBase.__has_mass_matrix(prob.f) && !(prob.f.mass_matrix isa LinearAlgebra.UniformScaling)
        OrdinaryDiffEqCore._initialize_dae!(integrator, prob, SciMLBase.CheckInit(), isinplace)
    end
end
