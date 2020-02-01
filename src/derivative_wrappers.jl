function derivative!(df::AbstractArray{<:Number}, f, x::Union{Number,AbstractArray{<:Number}}, fx::AbstractArray{<:Number}, integrator, grad_config)
    alg = unwrap_alg(integrator, true)
    tmp = length(x) # We calculate derivtive for all elements in gradient
    if alg_autodiff(alg)
        ForwardDiff.derivative!(df, f, fx, x, grad_config)
        integrator.destats.nf += 1
    else
        FiniteDiff.finite_difference_gradient!(df, f, x, grad_config, dir = diffdir(integrator))
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
      d = FiniteDiff.finite_difference_derivative(f, x, alg.diff_type, dir = diffdir(integrator))
      if alg.diff_type == Val{:central} || alg.diff_type == Val{:forward}
          tmp *= 2
      end
      integrator.destats.nf += tmp
      d
    end
end

jacobian_autodiff(f, x, odefun) = (ForwardDiff.derivative(f,x),1)
function jacobian_autodiff(f, x::AbstractArray, odefun)
  if DiffEqBase.has_colorvec(odefun)
    colorvec = odefun.colorvec
    sparsity = odefun.jac_prototype
    jac_prototype = nothing
  else
    colorvec = 1:length(x)
    sparsity = nothing
    jac_prototype = odefun.jac_prototype
  end
  maxcolor = maximum(colorvec)
  chunksize = getsize(default_chunk_size(maxcolor))
  num_of_chunks = Int(ceil(maxcolor / chunksize))
  (forwarddiff_color_jacobian(f,x,colorvec = colorvec, sparsity = sparsity,
   jac_prototype = jac_prototype),
   num_of_chunks)
end

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

jacobian_finitediff(f, x, diff_type, dir, colorvec, sparsity, jac_prototype) =
    (FiniteDiff.finite_difference_derivative(f, x, diff_type, eltype(x), dir = dir),2)
jacobian_finitediff(f, x::AbstractArray, diff_type, dir, colorvec, sparsity, jac_prototype) =
    (FiniteDiff.finite_difference_jacobian(f, x, diff_type, eltype(x), diff_type==Val{:forward} ? f(x) : similar(x),
      dir = dir, colorvec = colorvec, sparsity = sparsity, jac_prototype = jac_prototype),_nfcount(maximum(colorvec),diff_type))

function jacobian(f, x, integrator)
    alg = unwrap_alg(integrator, true)
    local tmp
    if alg_autodiff(alg)
      J, tmp = jacobian_autodiff(f, x, integrator.f)
    else
      if DiffEqBase.has_colorvec(integrator.f)
        colorvec = integrator.f.colorvec
        sparsity = integrator.f.jac_prototype
        jac_prototype = nothing
      else
        colorvec = 1:length(x)
        sparsity = nothing
        jac_prototype = integrator.f.jac_prototype
      end
      dir = diffdir(integrator)
      J, tmp = jacobian_finitediff(f, x, alg.diff_type, dir, colorvec, sparsity, jac_prototype)
    end
    integrator.destats.nf += tmp
    J
end

jacobian_finitediff_forward!(J,f,x,jac_config,forwardcache,integrator,colorvec)=
  (FiniteDiff.finite_difference_jacobian!(J,f,x,jac_config,forwardcache,
    dir=diffdir(integrator),colorvec=colorvec,sparsity=integrator.f.jac_prototype);maximum(colorvec))
jacobian_finitediff!(J,f,x,jac_config,integrator,colorvec)=
  (FiniteDiff.finite_difference_jacobian!(J,f,x,jac_config,
    dir=diffdir(integrator),colorvec=colorvec,sparsity=integrator.f.jac_prototype);2*maximum(colorvec))

function jacobian!(J::AbstractMatrix{<:Number}, f, x::AbstractArray{<:Number}, fx::AbstractArray{<:Number}, integrator::DiffEqBase.DEIntegrator, jac_config)
    alg = unwrap_alg(integrator, true)
    if alg_autodiff(alg)
      forwarddiff_color_jacobian!(J,f,x,jac_config)
      integrator.destats.nf += 1
    else
      colorvec = DiffEqBase.has_colorvec(integrator.f) ? integrator.f.colorvec : 1:length(x)
      isforward = alg.diff_type === Val{:forward}
      if isforward
        forwardcache = get_tmp_cache(integrator, alg, unwrap_cache(integrator, true))[2]
        f(forwardcache, x)
        integrator.destats.nf += 1
        tmp=jacobian_finitediff_forward!(J, f, x, jac_config, forwardcache, integrator, colorvec)
      else # not forward difference
        tmp=jacobian_finitediff!(J, f, x, jac_config, integrator, colorvec)
      end
      integrator.destats.nf += tmp
    end
    nothing
end

function DiffEqBase.build_jac_config(alg::Union{OrdinaryDiffEqAlgorithm,DAEAlgorithm},f,uf,du1,uprev,u,tmp,du2,::Val{transform}=Val(true)) where transform
  if !DiffEqBase.has_jac(f) && ((!transform && !DiffEqBase.has_Wfact(f)) || (transform && !DiffEqBase.has_Wfact_t(f)))
    if alg_autodiff(alg)
      if DiffEqBase.has_colorvec(f)
        colorvec = f.colorvec
        sparsity = f.jac_prototype
      else
        colorvec = 1:length(uprev)
        sparsity = nothing
      end
      jac_config = ForwardColorJacCache(uf,uprev,colorvec=colorvec,sparsity=sparsity)
    else
      if alg.diff_type != Val{:complex}
        jac_config = FiniteDiff.JacobianCache(tmp,du1,du2,alg.diff_type)
      else
        jac_config = FiniteDiff.JacobianCache(Complex{eltype(tmp)}.(tmp),Complex{eltype(du1)}.(du1),nothing,alg.diff_type,eltype(u))
      end
    end
  else
    jac_config = nothing
  end
  jac_config
end

get_chunksize(jac_config::ForwardDiff.JacobianConfig{T,V,N,D}) where {T,V,N,D} = N

function DiffEqBase.resize_jac_config!(jac_config::SparseDiffTools.ForwardColorJacCache, i)
  resize!(jac_config.fx, i)
  resize!(jac_config.dx, i)
  resize!(jac_config.t, i)
  ps = SparseDiffTools.adapt.(typeof(jac_config.dx),
                 SparseDiffTools.generate_chunked_partials(jac_config.dx,
                 1:length(jac_config.dx),Val(ForwardDiff.npartials(jac_config.t[1]))))
  resize!(jac_config.p, length(ps))
  jac_config.p .= ps
end

function DiffEqBase.resize_jac_config!(jac_config::FiniteDiff.JacobianCache, i)
  resize!(jac_config, i)
  jac_config
end

function resize_grad_config!(grad_config::ForwardDiff.DerivativeConfig, i)
  resize!(grad_config.duals, i)
  grad_config
end

function resize_grad_config!(grad_config::FiniteDiff.GradientCache, i)
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
      grad_config = FiniteDiff.GradientCache(du1,t,alg.diff_type)
    end
  else
    grad_config = nothing
  end
  grad_config
end
