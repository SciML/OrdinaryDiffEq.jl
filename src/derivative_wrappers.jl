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

function build_jac_config(alg,uf,du1,uprev,u,tmp,du2)
  if !has_jac(f)
    if alg_autodiff(alg)
      jac_config = ForwardDiff.JacobianConfig(uf,du1,uprev,ForwardDiff.Chunk{determine_chunksize(u,alg)}())
    else
      RealOrComplex = eltype(u) <: Complex ? Val{:Complex} : Val{:Real}
      if alg.diff_type != Val{:complex}
        jac_config = DiffEqDiffTools.JacobianCache(alg.diff_type,RealOrComplex,tmp,du1,du2)
      else
        jac_config = DiffEqDiffTools.JacobianCache(alg.diff_type,RealOrComplex,Complex{eltype(tmp)}.(tmp),Complex{eltype(du1)}.(du1),nothing)
      end
    end
  else
    jac_config = nothing
  end
  jac_config
end
