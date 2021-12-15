function derivative!(df::AbstractArray{<:Number}, f, x::Union{Number,AbstractArray{<:Number}}, fx::AbstractArray{<:Number}, integrator, grad_config)
    alg = unwrap_alg(integrator, true)
    tmp = length(x) # We calculate derivtive for all elements in gradient
    if alg_autodiff(alg)
        T = if standardtag(alg)
          typeof(ForwardDiff.Tag(OrdinaryDiffEqTag(),eltype(df)))
        else
          typeof(ForwardDiff.Tag(f,eltype(df)))
        end

        xdual = Dual{T,eltype(df),1}(convert(eltype(df),x),ForwardDiff.Partials((one(eltype(df)),)))
        f(grad_config,xdual)
        df .= first.(ForwardDiff.partials.(grad_config))
        integrator.destats.nf += 1
    else
        FiniteDiff.finite_difference_gradient!(df, f, x, grad_config, dir = diffdir(integrator))
        fdtype = alg_difftype(alg)
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
      d = FiniteDiff.finite_difference_derivative(f, x, alg_difftype(alg), dir = diffdir(integrator))
      if alg_difftype(alg) === Val{:central} || alg_difftype(alg) === Val{:forward}
          tmp *= 2
      end
      integrator.destats.nf += tmp
      d
    end
end

jacobian_autodiff(f, x, odefun, alg) = (ForwardDiff.derivative(f,x),1, alg)
function jacobian_autodiff(f, x::AbstractArray, odefun, alg)
  jac_prototype = odefun.jac_prototype
  sparsity,colorvec = sparsity_colorvec(odefun,x)
  maxcolor = maximum(colorvec)
  chunk_size = get_chunksize(alg)===Val(0) ? nothing : get_chunksize(alg)
  num_of_chunks = chunk_size===nothing ? Int(ceil(maxcolor / getsize(ForwardDiff.pickchunksize(maxcolor)))) :
                                        Int(ceil(maxcolor / _unwrap_val(chunk_size)))
  (forwarddiff_color_jacobian(f,x,colorvec = colorvec, sparsity = sparsity,
                              jac_prototype = jac_prototype, chunksize=chunk_size),
   num_of_chunks)
end

function _nfcount(N,::Type{diff_type}) where diff_type
  if diff_type===Val{:complex}
    tmp = N
  elseif diff_type===Val{:forward}
    tmp = N + 1
  else
    tmp = 2N
  end
  tmp
end

jacobian_finitediff(f, x, ::Type{diff_type}, dir, colorvec, sparsity, jac_prototype) where {diff_type} =
    (FiniteDiff.finite_difference_derivative(f, x, diff_type, eltype(x), dir = dir),2)
function jacobian_finitediff(f, x::AbstractArray, ::Type{diff_type}, dir, colorvec, sparsity, jac_prototype) where {diff_type}
  f_in = diff_type === Val{:forward} ? f(x) : similar(x)
  ret_eltype = eltype(f_in)
  J = FiniteDiff.finite_difference_jacobian(f, x, diff_type, ret_eltype, f_in,
                                            dir = dir, colorvec = colorvec, sparsity = sparsity, jac_prototype = jac_prototype)
  return J, _nfcount(maximum(colorvec),diff_type)
end
function jacobian(f, x, integrator)
    alg = unwrap_alg(integrator, true)
    local tmp
    if alg_autodiff(alg)
      J, tmp = jacobian_autodiff(f, x, integrator.f, alg)
    else
      jac_prototype = integrator.f.jac_prototype
      sparsity,colorvec = sparsity_colorvec(integrator.f,x)
      dir = diffdir(integrator)
      J, tmp = jacobian_finitediff(f, x, alg_difftype(alg), dir, colorvec, sparsity, jac_prototype)
    end
    integrator.destats.nf += tmp
    J
end

jacobian_finitediff_forward!(J,f,x,jac_config,forwardcache,integrator)=
  (FiniteDiff.finite_difference_jacobian!(J,f,x,jac_config,forwardcache,
    dir=diffdir(integrator));maximum(jac_config.colorvec))
jacobian_finitediff!(J,f,x,jac_config,integrator)=
  (FiniteDiff.finite_difference_jacobian!(J,f,x,jac_config,
    dir=diffdir(integrator));2*maximum(jac_config.colorvec))

function jacobian!(J::AbstractMatrix{<:Number}, f, x::AbstractArray{<:Number}, fx::AbstractArray{<:Number}, integrator::DiffEqBase.DEIntegrator, jac_config)
    alg = unwrap_alg(integrator, true)
    if alg_autodiff(alg)
      forwarddiff_color_jacobian!(J,f,x,jac_config)
      integrator.destats.nf += 1
    else
      isforward = alg_difftype(alg) === Val{:forward}
      if isforward
        forwardcache = get_tmp_cache(integrator, alg, unwrap_cache(integrator, true))[2]
        f(forwardcache, x)
        integrator.destats.nf += 1
        tmp=jacobian_finitediff_forward!(J, f, x, jac_config, forwardcache, integrator)
      else # not forward difference
        tmp=jacobian_finitediff!(J, f, x, jac_config, integrator)
      end
      integrator.destats.nf += tmp
    end
    nothing
end

function DiffEqBase.build_jac_config(alg,f::F1,uf::F2,du1,uprev,u,tmp,du2,::Val{transform}=Val(true)) where {transform,F1,F2}
  if !DiffEqBase.has_jac(f) && ((!transform && !DiffEqBase.has_Wfact(f)) || (transform && !DiffEqBase.has_Wfact_t(f)))
    jac_prototype = f.jac_prototype
    sparsity,colorvec = sparsity_colorvec(f,u)
    if alg_autodiff(alg)
      _chunksize = get_chunksize(alg)===Val(0) ? nothing : get_chunksize(alg) # SparseDiffEq uses different convection...

      T = if standardtag(alg)
        typeof(ForwardDiff.Tag(OrdinaryDiffEqTag(),eltype(u)))
      else
        typeof(ForwardDiff.Tag(uf,eltype(u)))
      end

      jac_config = ForwardColorJacCache(uf,uprev,_chunksize;colorvec=colorvec,sparsity=sparsity,tag=T)
    else
      if alg_difftype(alg) !== Val{:complex}
        jac_config = FiniteDiff.JacobianCache(tmp,du1,du2,alg_difftype(alg),colorvec=colorvec,sparsity=sparsity)
      else
        jac_config = FiniteDiff.JacobianCache(Complex{eltype(tmp)}.(tmp),Complex{eltype(du1)}.(du1),nothing,alg_difftype(alg),eltype(u),colorvec=colorvec,sparsity=sparsity)
      end
    end
  else
    jac_config = nothing
  end
  jac_config
end

get_chunksize(jac_config::ForwardDiff.JacobianConfig{T,V,N,D}) where {T,V,N,D} = Val(N) # don't degrade compile time information to runtime information

function DiffEqBase.resize_jac_config!(jac_config::SparseDiffTools.ForwardColorJacCache, i)
  resize!(jac_config.fx, i)
  resize!(jac_config.dx, i)
  resize!(jac_config.t, i)
  ps = SparseDiffTools.adapt.(DiffEqBase.parameterless_type(jac_config.dx),
                 SparseDiffTools.generate_chunked_partials(jac_config.dx,
                 1:length(jac_config.dx),Val(ForwardDiff.npartials(jac_config.t[1]))))
  resize!(jac_config.p, length(ps))
  jac_config.p .= ps
end

function DiffEqBase.resize_jac_config!(jac_config::FiniteDiff.JacobianCache, i)
  resize!(jac_config, i)
  jac_config
end

function resize_grad_config!(grad_config::AbstractArray, i)
  resize!(grad_config, i)
  grad_config
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

function build_grad_config(alg,f::F1,tf::F2,du1,t) where {F1,F2}
  if !DiffEqBase.has_tgrad(f)
    if alg_autodiff(alg)

      T = if standardtag(alg)
        typeof(ForwardDiff.Tag(OrdinaryDiffEqTag(),eltype(du1)))
      else
        typeof(ForwardDiff.Tag(f,eltype(du1)))
      end

      if du1 isa Array
        dualt = Dual{T,eltype(du1),1}(first(du1)*t, ForwardDiff.Partials((one(eltype(du1)),)))
        grad_config = similar(du1,typeof(dualt))
        fill!(grad_config,false)
      else
        grad_config = Dual{T,eltype(du1),1}.(du1, (ForwardDiff.Partials((one(eltype(du1)),)),)) .* false
      end
    else
      grad_config = FiniteDiff.GradientCache(du1,t,alg_difftype(alg))
    end
  else
    grad_config = nothing
  end
  grad_config
end

function sparsity_colorvec(f,x)
  sparsity = f.sparsity
  colorvec = DiffEqBase.has_colorvec(f) ? f.colorvec :
              (isnothing(sparsity) ? (1:length(x)) : matrix_colors(sparsity))
  sparsity,colorvec
end
