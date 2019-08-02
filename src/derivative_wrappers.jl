function derivative!(df::AbstractArray{<:Number}, f, x::Union{Number,AbstractArray{<:Number}}, fx::AbstractArray{<:Number}, integrator, grad_config)
    alg = unwrap_alg(integrator, true)
    tmp = length(x) # We calculate derivtive for all elements in gradient
    if alg_autodiff(alg)
        ForwardDiff.derivative!(df, f, fx, x, grad_config)
        integrator.destats.nf += 1
    else
        DiffEqDiffTools.finite_difference_gradient!(df, f, x, grad_config, dir = diffdir(integrator))
        fdtype = alg.diff_type
        if fdtype == Val{:forward} || fdtype == Val{:central}
            tmp *= 2
            if eltype(df)<:Complex
              tmp *= 2
            end
        end
        integrator.destats.nf += tmp
    end
    nothing
end

function derivative(f, x::Union{Number,AbstractArray{<:Number}},
                    integrator)
    local d
    tmp = length(x) # We calculate derivtive for all elements in gradient
    alg = unwrap_alg(integrator, true)
    if alg_autodiff(alg)
      integrator.destats.nf += 1
      d = ForwardDiff.derivative(f, x)
    else
      d = DiffEqDiffTools.finite_difference_derivative(f, x, alg.diff_type, dir = diffdir(integrator))
      if alg.diff_type == Val{:central} || alg.diff_type == Val{:forward}
          tmp *= 2
      end
      integrator.destats.nf += tmp
      d
    end
end

function jacobian(f, x, integrator)
    alg = unwrap_alg(integrator, true)
    local tmp
    if alg_autodiff(alg)
      if DiffEqBase.has_colorvec(integrator.f)
        J,tmp = jacobian_autodiff(f, x, integrator)
      else
        J,tmp = jacobian_autodiff(f, x)
      end
    else
      if DiffEqBase.has_colorvec(integrator.f)
        J,tmp = jacobian_finitediff(f, x, alg.diff_type, integrator, integrator.f.colorvec)
      else
        J,tmp = jacobian_finitediff(f, x, alg.diff_type, integrator)
      end
    end
    integrator.destats.nf += tmp
    J
end

jacobian_autodiff(f, x) = (ForwardDiff.derivative(f,x),1)
jacobian_autodiff(f, x::AbstractArray) = (ForwardDiff.jacobian(f, x),1)
function jacobian_autodiff(f, x::AbstractArray, integrator)
  colorvec=integrator.f.colorvec
  jac=integrator.f.jac_prototype
  J=jac isa SparseMatrixCSC ? similar(jac) : fill(0.,size(jac))
  (forwarddiff_color_jacobian!(J,f,x,colorvec=colorvec,sparsity=jac),1)
end
#jacobian_autodiff(f, x::AbstractArray, colorvec) = (ForwardDiff.jacobian(f, x, colorvec = colorvec),1)

function _nfcount(N,diff_type)
  if diff_type==Val{:complex}
    tmp = N
  elseif diff_type==Val{:forward}
    tmp = N + 1
  else
    tmp = 2N
  end
  tmp
end

jacobian_finitediff(f, x, diff_type, integrator) =
    (DiffEqDiffTools.finite_difference_derivative(f, x, diff_type, eltype(x), dir = diffdir(integrator)),2)
jacobian_finitediff(f, x::AbstractArray, diff_type, integrator) =
    (DiffEqDiffTools.finite_difference_jacobian(f, x, diff_type, eltype(x), Val{false}, 
      dir = diffdir(integrator)),_nfcount(length(x),diff_type))
jacobian_finitediff(f, x::AbstractArray, diff_type, integrator, colorvec) =
    (DiffEqDiffTools.finite_difference_jacobian(f, x, diff_type, eltype(x), Val{false}, 
      dir = diffdir(integrator), colorvec = colorvec, sparsity = integrator.f.jac_prototype),_nfcount(maximum(colorvec),diff_type))

jacobian_finitediff_forward!(J,f,x,jac_config,forwardcache,integrator)=(DiffEqDiffTools.finite_difference_jacobian!(J,f,x,jac_config,forwardcache,dir=diffdir(integrator));length(x))
jacobian_finitediff_forward!(J,f,x,jac_config,forwardcache,integrator,colorvec)=
  (DiffEqDiffTools.finite_difference_jacobian!(J,f,x,jac_config,forwardcache,
    dir=diffdir(integrator),colorvec=colorvec,sparsity=integrator.f.jac_prototype);maximum(colorvec))
jacobian_finitediff!(J,f,x,jac_config,integrator)=(DiffEqDiffTools.finite_difference_jacobian!(J,f,x,jac_config,dir=diffdir(integrator));2*length(x))
jacobian_finitediff!(J,f,x,jac_config,integrator,colorvec)=
  (DiffEqDiffTools.finite_difference_jacobian!(J,f,x,jac_config,
    dir=diffdir(integrator),colorvec=colorvec,sparsity=integrator.f.jac_prototype);2*maximum(colorvec))
jacobian_autodiff!(J,f,x,jac_config)=forwarddiff_color_jacobian!(J,f,x,jac_config)#J::SparseMatrixCSC
jacobian_autodiff!(J,f,fx,x,jac_config)=ForwardDiff.jacobian!(J,f,fx,x,jac_config)

function jacobian!(J::AbstractMatrix{<:Number}, f, x::AbstractArray{<:Number}, fx::AbstractArray{<:Number}, integrator::DiffEqBase.DEIntegrator, jac_config)
    alg = unwrap_alg(integrator, true)
    if alg_autodiff(alg)
      if DiffEqBase.has_colorvec(integrator.f)
        jacobian_autodiff!(J,f,x,jac_config)
      else
        jacobian_autodiff!(J,f,fx,x,jac_config)
      end
      integrator.destats.nf += 1
    else
      isforward = alg.diff_type === Val{:forward}
      if isforward
        forwardcache = get_tmp_cache(integrator, alg, unwrap_cache(integrator, true))[2]
        f(forwardcache, x)
        integrator.destats.nf += 1
        if DiffEqBase.has_colorvec(integrator.f)
          tmp=jacobian_finitediff_forward!(J, f, x, jac_config, forwardcache, integrator, integrator.f.colorvec)
        else
          tmp=jacobian_finitediff_forward!(J, f, x, jac_config, forwardcache, integrator)
        end
      else # not forward difference
        if DiffEqBase.has_colorvec(integrator.f)
          tmp=jacobian_finitediff!(J, f, x, jac_config, integrator, integrator.f.colorvec)
        else
          tmp=jacobian_finitediff!(J, f, x, jac_config, integrator)
        end
      end
      integrator.destats.nf += tmp
    end
    nothing
end

function DiffEqBase.build_jac_config(alg::OrdinaryDiffEqAlgorithm,f,uf,du1,uprev,u,tmp,du2)
  if !DiffEqBase.has_jac(f)
    if alg_autodiff(alg)
      if DiffEqBase.has_colorvec(f)
        jac_config = ForwardColorJacCache(uf,uprev,colorvec=f.colorvec,sparsity=f.jac_prototype)
      else
        jac_config = ForwardDiff.JacobianConfig(uf,du1,uprev,ForwardDiff.Chunk{determine_chunksize(u,alg)}())
      end
    else
      if alg.diff_type != Val{:complex}
        jac_config = DiffEqDiffTools.JacobianCache(tmp,du1,du2,alg.diff_type)
      else
        jac_config = DiffEqDiffTools.JacobianCache(Complex{eltype(tmp)}.(tmp),Complex{eltype(du1)}.(du1),nothing,alg.diff_type,eltype(u))
      end
    end
  else
    jac_config = nothing
  end
  jac_config
end

get_chunksize(jac_config::ForwardDiff.JacobianConfig{T,V,N,D}) where {T,V,N,D} = N

function DiffEqBase.resize_jac_config!(jac_config::ForwardDiff.JacobianConfig, i)
  for j in eachindex(jac_config.duals)
    resize!(jac_config.duals[j], i)
  end
  jac_config
end

function DiffEqBase.resize_jac_config!(jac_config::DiffEqDiffTools.JacobianCache, i)
  resize!(jac_config, i)
  jac_config
end

function resize_grad_config!(grad_config::ForwardDiff.DerivativeConfig, i)
  resize!(grad_config.duals, i)
  grad_config
end

function resize_grad_config!(grad_config::DiffEqDiffTools.GradientCache, i)
  @unpack fx, c1, c2 = grad_config
  fx !== nothing && resize!(fx, i)
  c1 !== nothing && resize!(c1, i)
  c2 !== nothing && resize!(c2, i)
  grad_config
end

function build_grad_config(alg,f,tf,du1,t)
  if !DiffEqBase.has_tgrad(f)
    if alg_autodiff(alg)
      grad_config = ForwardDiff.DerivativeConfig(tf,du1,t)
    else
      grad_config = DiffEqDiffTools.GradientCache(du1,t,alg.diff_type)
    end
  else
    grad_config = nothing
  end
  grad_config
end
