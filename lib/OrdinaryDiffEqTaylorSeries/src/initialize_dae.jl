# Initialize DAE for DAETS algorithm
function OrdinaryDiffEqCore._initialize_dae!(integrator::OrdinaryDiffEqCore.ODEIntegrator{<:DAETS}, 
                                           prob::DAEProblem, 
                                           alg::OrdinaryDiffEqCore.DefaultInit, 
                                           ::Val{true})
    @unpack f, u0, du0, p, tspan = prob
    t0 = first(tspan) 
    
    #integrator
    integrator.u = u0
    integrator.uprev = copy(u0)
    integrator.t = t0
    
    #cache
    initialize!(integrator, integrator.cache)
    integrator.k[1] .= du0
    
    # Step: run structural analysis or find jacobian or something TODO: fix this
    # Σ, modified_eqs = structural_analysis(f, u0, tspan, p)
    return nothing
end 