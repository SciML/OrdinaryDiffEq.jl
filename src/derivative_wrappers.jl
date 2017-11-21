function derivative!(df::AbstractArray{<:Number}, f, x::Union{Number,AbstractArray{<:Number}}, fx::AbstractArray{<:Number}, integrator::DEIntegrator)
    if alg_autodiff(integrator.alg)
        ForwardDiff.derivative!(df, f, fx, x)
    else
        RealOrComplex = eltype(integrator.u) <: Complex ? Val{:Complex} : Val{:Real}
        DiffEqDiffTools.finite_difference!(df, f, x, integrator.alg.diff_type, RealOrComplex, fx)
    end
    nothing
end

function jacobian!(J::AbstractMatrix{<:Number}, f, x::AbstractArray{<:Number}, fx::AbstractArray{<:Number}, integrator::DEIntegrator, jac_config)
    if alg_autodiff(integrator.alg)
      ForwardDiff.jacobian!(J, f, fx, x, jac_config)
    else
      DiffEqDiffTools.finite_difference_jacobian!(J, f, x, jac_config)
    end
    nothing
end
