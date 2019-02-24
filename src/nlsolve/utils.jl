DiffEqBase.@def iipnlsolve begin
  # unpack cache of non-linear solver
  _nlcache = alg.nlsolve.cache
  @unpack max_iter,min_iter = _nlcache

  # define additional fields of cache of non-linear solver
  z = similar(u); dz = similar(u); tmp = similar(u); b = similar(u)
  k = zero(rate_prototype)

  # define tolerances
  uToltype = real(uBottomEltypeNoUnits)
  κ = _nlcache.κ === nothing ? uToltype(1//100) : uToltype(_nlcache.κ)
  tol = _nlcache.tol === nothing ? uToltype(min(0.03,first(reltol)^(0.5))) : uToltype(_nlcache.tol)
  ηold = one(uToltype)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    nf = nlsolve_f(f, alg)

    # check if `nf` is linear
    islin = f isa Union{ODEFunction,SplitFunction} && islinear(nf.f)

    if islin
      # get the operator
      J = nf.f
      W = WOperator(f.mass_matrix, dt, J, true)
    else
      if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype !== nothing
        W = WOperator(f, dt, true)
        J = nothing # is J = W.J better?
      else
        J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
        W = similar(J)
      end
    end

    nlcache = NLNewtonCache(κ,tol,min_iter,max_iter,10000,_nlcache.new_W,z,W,γ,c,ηold,dz,tmp,b,k)
  elseif alg.nlsolve isa NLFunctional
    z₊ = similar(z)

    nlcache = NLFunctionalCache(κ,tol,min_iter,max_iter,10000,z,γ,c,ηold,z₊,dz,tmp,b,k)
  elseif alg.nlsolve isa NLAnderson
    z₊ = similar(z)
    zs = [zero(vec(z)) for i in 1:alg.nlsolve.n+1]
    gs = [zero(vec(z)) for i in 1:alg.nlsolve.n+1]
    alphas = Array{eltype(z)}(undef, alg.nlsolve.n + 1)
    residuals = Array{eltype(z)}(undef, length(z), alg.nlsolve.n + 1)

    nlcache = NLAndersonCache(κ,tol,min_iter,max_iter,10000,z,γ,c,ηold,alphas,residuals,z₊,dz,tmp,b,k,zs,gs)
  end

  # create non-linear solver
  _nlsolve = typeof(alg.nlsolve).name.wrapper
  nlsolve = _nlsolve{true, typeof(nlcache)}(nlcache)

  # define additional fields of cache
  fsalfirst = zero(rate_prototype)
  if alg.nlsolve isa NLNewton
    if islin
      du1 = rate_prototype
      uf = nothing
      jac_config = nothing
      linsolve = alg.linsolve(Val{:init},nf,u)
    else
      du1 = zero(rate_prototype)
      # if the algorithm specializes on split problems the use `nf`
      uf = DiffEqDiffTools.UJacobianWrapper(nf,t,p)
      jac_config = build_jac_config(alg,nf,uf,du1,uprev,u,tmp,dz)
      linsolve = alg.linsolve(Val{:init},uf,u)
    end
  else
    J = nothing
    W = nothing
    du1 = rate_prototype
    uf = nothing
    jac_config = nothing
    linsolve = nothing
  end
end

DiffEqBase.@def oopnlsolve begin
  # unpack cache of non-linear solver
  _nlcache = alg.nlsolve.cache
  @unpack max_iter,min_iter = _nlcache

  # define additional fields of cache of non-linear solver (all aliased)
  z = uprev; dz = z; tmp = z; b = z; k = rate_prototype

  # define tolerances
  uToltype = real(uBottomEltypeNoUnits)
  κ = _nlcache.κ === nothing ? uToltype(1//100) : uToltype(_nlcache.κ)
  tol = _nlcache.tol === nothing ? uToltype(min(0.03,first(reltol)^(0.5))) : uToltype(_nlcache.tol)
  ηold = one(uToltype)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    nf = nlsolve_f(f, alg)

    # only use `nf` if the algorithm specializes on split eqs
    uf = DiffEqDiffTools.UDerivativeWrapper(nf,t,p)

    islin = f isa Union{ODEFunction,SplitFunction} && islinear(nf.f)
    if islin || DiffEqBase.has_jac(f)
      # get the operator
      J = islin ? nf.f : f.jac(uprev, p, t)
      if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
        J = DiffEqArrayOperator(J)
      end
      W = WOperator(f.mass_matrix, dt, J, false)
    else
      #if DiffEqBase.has_jac(f)
      #  J = f.jac(uprev, p, t)
      #else
      #  if alg_autodiff(alg)
      #    J = jacobian_autodiff(uf, uprev)
      #  else
      #    J = jacobian_finitediff(uf, uprev, alg.diff_type)
      #  end
      #end
      #W = J isa Number ? J : lu(J; check=false)
      W = u isa Number ? u : LU{LinearAlgebra.lutype(eltype(u))}(Matrix{uEltypeNoUnits}(undef, 0, 0),
                                                                 zeros(LinearAlgebra.BlasInt, min(size(J)...)),
                                                                 0)
    end

    nlcache = NLNewtonCache(κ,tol,min_iter,max_iter,10000,_nlcache.new_W,z,W,γ,c,ηold,dz,tmp,b,k)
  elseif alg.nlsolve isa NLFunctional
    uf = nothing

    nlcache = NLFunctionalCache(κ,tol,min_iter,max_iter,10000,z,γ,c,ηold,z,dz,tmp,b,k)
  elseif alg.nlsolve isa NLAnderson
    uf = nothing
    zs = fill(_vec(z), alg.nlsolve.n + 1)
    gs = fill(_vec(z), alg.nlsolve.n + 1)
    alphas = nothing
    residuals = Array{eltype(z)}(undef, length(z), alg.nlsolve.n + 1)

    nlcache = NLAndersonCache(κ,tol,min_iter,max_iter,10000,z,γ,c,ηold,alphas,residuals,z,dz,tmp,b,k,zs,gs)
  end

  # create non-linear solver
  _nlsolve = typeof(alg.nlsolve).name.wrapper
  nlsolve = _nlsolve{false, typeof(nlcache)}(nlcache)
end
