get_status(nlsolver::NLSolver) = nlsolver.status
nlsolvefail(nlsolver::NLSolver) = nlsolvefail(get_status(nlsolver))
nlsolvefail(status::NLStatus) = Int8(status) < 0

isnewton(nlsolver::NLSolver) = isnewton(nlsolver.cache)
isnewton(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}) = true
isnewton(nlcache::AbstractNLSolverCache) = false

set_new_W!(nlsolver::NLSolver, val::Bool)::Bool = set_new_W!(nlsolver.cache, val)
set_new_W!(nlcache::NLNewtonCache, val::Bool)::Bool = nlcache.new_W = val
set_new_W!(nlcache::AbstractNLSolverCache, val::Bool)::Bool = val

get_W(nlsolver::NLSolver) = get_W(nlsolver.cache)
get_W(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}) = nlcache.W

set_W!(nlsolver::NLSolver, W) = set_W!(nlsolver.cache, W)
set_W!(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}, W) = (nlcache.W = W; W)

set_W_dt!(nlsolver::NLSolver, W_dt) = set_W_dt!(nlsolver.cache, W_dt)
set_W_dt!(nlcache::NLNewtonCache, W_dt) = (nlcache.W_dt = W_dt; W_dt)
set_W_dt!(nlcache::NLNewtonConstantCache, W_dt) = W_dt

DiffEqBase.@def iipnlsolve begin
  @unpack κ, fast_convergence_cutoff = alg.nlsolve

  # define additional fields of cache of non-linear solver
  z = similar(u); dz = similar(u); tmp = similar(u); b = similar(u)
  k = zero(rate_prototype); ustep = similar(u)
  tstep = zero(t)
  atmp = similar(u, uEltypeNoUnits)

  uTolType = real(uBottomEltypeNoUnits)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    nf = nlsolve_f(f, alg)

    if islinear(f)
      # get the operator
      J = nf.f
      W = WOperator(f.mass_matrix, dt, J, true)
    else
      if DiffEqBase.has_jac(f) && !DiffEqBase.has_Wfact(f) && f.jac_prototype !== nothing
        W = WOperator(f, dt, true)
        J = nothing # is J = W.J better?
      else
        J = false .* vec(u) .* vec(u)'
        W = similar(J)
      end
    end

    tType = typeof(t)
    invγdt = inv(oneunit(t) * one(uTolType))

    nlcache = NLNewtonCache(ustep,tstep,atmp,true,W,J,tType(dt),invγdt,tType(alg.nlsolve.new_W_dt_cutoff))
  elseif alg.nlsolve isa NLFunctional
    nlcache = NLFunctionalCache(ustep,tstep,atmp)
  elseif alg.nlsolve isa NLAnderson
    max_history = min(alg.nlsolve.max_history, alg.nlsolve.max_iter, length(z))
    Δz₊s = [zero(z) for i in 1:max_history]
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)
    dzold = zero(z)
    z₊old = zero(z)

    nlcache = NLAndersonCache(ustep,tstep,atmp,dzold,z₊old,Δz₊s,Q,R,γs,0,alg.nlsolve.aa_start,alg.nlsolve.droptol)
  end

  # define additional fields of cache
  fsalfirst = zero(rate_prototype)
  if alg.nlsolve isa NLNewton
    if islinear(f)
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
    # TODO: check if the solver is iterative
    weight = similar(u)
  else
    J = nothing
    W = nothing
    du1 = rate_prototype
    uf = nothing
    jac_config = nothing
    linsolve = nothing
    weight = z
  end

  # create non-linear solver
  nlsolver = NLSolver{true,typeof(z),typeof(k),uTolType,typeof(κ),typeof(γ),typeof(c),typeof(du1),typeof(uf),typeof(jac_config),typeof(linsolve),typeof(fast_convergence_cutoff),typeof(nlcache)}(z,dz,tmp,b,k,one(uTolType),κ,γ,c,alg.nlsolve.max_iter,10000,Convergence,fast_convergence_cutoff,du1,uf,jac_config,linsolve,weight,nlcache)
end

DiffEqBase.@def oopnlsolve begin
  @unpack κ, fast_convergence_cutoff = alg.nlsolve

  # define additional fields of cache of non-linear solver (all aliased)
  z = uprev; dz = z; tmp = z; b = z; k = rate_prototype
  tstep = zero(t)

  uTolType = real(uBottomEltypeNoUnits)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    nf = nlsolve_f(f, alg)

    # only use `nf` if the algorithm specializes on split eqs
    uf = DiffEqDiffTools.UDerivativeWrapper(nf,t,p)

    if islinear(f) || DiffEqBase.has_jac(f)
      # get the operator
      J = islinear(f) ? nf.f : f.jac(uprev, p, t)
      if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
        J = DiffEqArrayOperator(J)
      end
      W = WOperator(f.mass_matrix, dt, J, false)
    else
      # https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/pull/672
      if u isa StaticArray
        # get a "fake" `J`
        J = if u isa AbstractMatrix && size(u, 1) > 1 # `u` is already a matrix
          u
        elseif size(u, 1) == 1 # `u` is a row vector
          vcat(u, u)
        else # `u` is a column vector
          hcat(u, u)
        end
        W = lu(J)
      else
        W = u isa Number ? u : LU{LinearAlgebra.lutype(uEltypeNoUnits)}(Matrix{uEltypeNoUnits}(undef, 0, 0),
                                                                        Vector{LinearAlgebra.BlasInt}(undef, 0),
                                                                        zero(LinearAlgebra.BlasInt))
      end
    end

    invγdt = inv(oneunit(t) * one(uTolType))

    nlcache = NLNewtonConstantCache(tstep,W,J,invγdt,typeof(t)(alg.nlsolve.new_W_dt_cutoff))
  elseif alg.nlsolve isa NLFunctional
    uf = nothing

    nlcache = NLFunctionalConstantCache(tstep)
  elseif alg.nlsolve isa NLAnderson
    uf = nothing

    max_history = min(alg.nlsolve.max_history, alg.nlsolve.max_iter, length(z))
    Δz₊s = Vector{typeof(z)}(undef, max_history)
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)
    dzold = u
    z₊old = u

    nlcache = NLAndersonConstantCache(tstep,dzold,z₊old,Δz₊s,Q,R,γs,0,alg.nlsolve.aa_start,alg.nlsolve.droptol)
  end

  # create non-linear solver
  nlsolver = NLSolver{false,typeof(z),typeof(k),uTolType,typeof(κ),typeof(γ),typeof(c),Nothing,typeof(uf),Nothing,Nothing,typeof(fast_convergence_cutoff),typeof(nlcache)}(z,dz,tmp,b,k,one(uTolType),κ,γ,c,alg.nlsolve.max_iter,10000,Convergence,fast_convergence_cutoff,nothing,uf,nothing,nothing,z,nlcache)
end

# No J version
function iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  iipnlsolve(alg,u,uprev,p,t,dt,f,W,nothing,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
end

function iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @unpack κ, fast_convergence_cutoff = alg.nlsolve

  # define additional fields of cache of non-linear solver
  z = similar(u); dz = similar(u); tmp = similar(u); b = similar(u)
  k = zero(rate_prototype); ustep = similar(u)
  tstep = zero(t)
  atmp = similar(u, uEltypeNoUnits)

  uTolType = real(uBottomEltypeNoUnits)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    tType = typeof(t)
    invγdt = inv(oneunit(t) * one(uTolType))

    nlcache = NLNewtonCache(ustep,tstep,atmp,true,W,J,tType(dt),invγdt,tType(alg.nlsolve.new_W_dt_cutoff))
  elseif alg.nlsolve isa NLFunctional
    nlcache = NLFunctionalCache(ustep,tstep,atmp)
  elseif alg.nlsolve isa NLAnderson
    max_history = min(alg.nlsolve.max_history, alg.nlsolve.max_iter, length(z))
    Δz₊s = [zero(z) for i in 1:max_history]
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)
    dzold = zero(z)
    z₊old = zero(z)

    nlcache = NLAndersonCache(ustep,tstep,atmp,dzold,z₊old,Δz₊s,Q,R,γs,0,alg.nlsolve.aa_start,alg.nlsolve.droptol)
  end

  # define additional fields of cache
  if alg.nlsolve isa NLNewton
    nf = nlsolve_f(f, alg)

    if islinear(f)
      du1 = rate_prototype
      uf = nothing
      jac_config = nothing
      linsolve = alg.linsolve(Val{:init},nf,u)
    else
      du1 = zero(rate_prototype)
      # if the algorithm specializes on split problems the use `nf`
      # we pass this `alg` here just for identification purpose, because get_uf would be overloaded in different repos
      uf = iip_get_uf(alg,nf,t,p)
      jac_config = build_jac_config(alg,nf,uf,du1,uprev,u,tmp,dz)
      linsolve = alg.linsolve(Val{:init},uf,u)
    end
    # TODO: check if the solver is iterative
    weight = similar(u)
  else
    du1 = rate_prototype
    uf = nothing
    jac_config = nothing
    linsolve = nothing
    weight = z
  end

  # create non-linear solver
  nlsolver = NLSolver{true,typeof(z),typeof(k),uTolType,typeof(κ),typeof(γ),typeof(c),typeof(du1),typeof(uf),typeof(jac_config),typeof(linsolve),typeof(fast_convergence_cutoff),typeof(nlcache)}(z,dz,tmp,b,k,one(uTolType),κ,γ,c,alg.nlsolve.max_iter,10000,Convergence,fast_convergence_cutoff,du1,uf,jac_config,linsolve,weight,nlcache)
end

DiffEqBase.@def getiipnlsolvefields begin
  @unpack z,dz,tmp,k,uf,du1,jac_config,linsolve = nlsolver
  b = nlsolver.ztmp
  fsalfirst = zero(rate_prototype)
end

# No J version
function oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  oopnlsolve(alg,u,uprev,p,t,dt,f,W,nothing,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)  
end

function oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @unpack κ, fast_convergence_cutoff = alg.nlsolve

  # define additional fields of cache of non-linear solver (all aliased)
  z = uprev; dz = z; tmp = z; b = z; k = rate_prototype
  tstep = zero(t)

  uTolType = real(uBottomEltypeNoUnits)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    nf = nlsolve_f(f, alg)
    # only use `nf` if the algorithm specializes on split eqs
    uf = oop_get_uf(alg,nf,t,p)

    invγdt = inv(oneunit(t) * one(uTolType))

    nlcache = NLNewtonConstantCache(tstep,W,J,invγdt,typeof(t)(alg.nlsolve.new_W_dt_cutoff))
  elseif alg.nlsolve isa NLFunctional
    uf = nothing

    nlcache = NLFunctionalConstantCache(tstep)
  elseif alg.nlsolve isa NLAnderson
    uf = nothing

    max_history = min(alg.nlsolve.max_history, alg.nlsolve.max_iter, length(z))
    Δz₊s = Vector{typeof(z)}(undef, max_history)
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)
    dzold = u
    z₊old = u

    nlcache = NLAndersonConstantCache(tstep,dzold,z₊old,Δz₊s,Q,R,γs,0,alg.nlsolve.aa_start,alg.nlsolve.droptol)
  end

  # create non-linear solver
  nlsolver = NLSolver{false,typeof(z),typeof(k),uTolType,typeof(κ),typeof(γ),typeof(c),Nothing,typeof(uf),Nothing,Nothing,typeof(fast_convergence_cutoff),typeof(nlcache)}(z,dz,tmp,b,k,one(uTolType),κ,γ,c,alg.nlsolve.max_iter,10000,Convergence,fast_convergence_cutoff,nothing,uf,nothing,nothing,z,nlcache)
end

DiffEqBase.@def getoopnlsolvefields begin
  uf = nlsolver.uf
end

function nlsolve_resize!(integrator::DiffEqBase.DEIntegrator, i::Int)
  if !isdefined(integrator.cache, :nlsolver)
    return nothing
  end
  alg = integrator.alg; nlsolver = integrator.cache.nlsolver
  if nlsolver isa AbstractArray
    for idx in eachindex(nlsolver) # looping because we may have multiple nlsolver for threaded case
      _nlsolver = nlsolver[idx]
      @unpack z,dz,tmp,ztmp,k,du1,uf,jac_config,linsolve,weight,cache = _nlsolver
      # doubt: if these fields are always going to be in alg cache too, then we shouldnt do this here.
      # double resize doesn't do any bad I think though
      resize!(z,i)
      resize!(dz,i)
      resize!(tmp,i)
      resize!(ztmp,i)
      resize!(k,i)
      resize!(du1,i)
      if jac_config !== nothing
        _nlsolver.jac_config = resize_jac_config!(jac_config, i)
      end
      resize!(weight, i)
      nlsolve_cache_resize!(cache,alg,i)
    end
  else
    @unpack z,dz,tmp,ztmp,k,du1,uf,jac_config,linsolve,weight,cache = nlsolver
    resize!(z,i)
    resize!(dz,i)
    resize!(tmp,i)
    resize!(ztmp,i)
    resize!(k,i)
    resize!(du1,i)
    if jac_config !== nothing
      nlsolver.jac_config = resize_jac_config!(jac_config,i)
    end
    resize!(weight, i)
    nlsolve_cache_resize!(cache,alg,i)
  end
  nothing
end

function nlsolve_cache_resize!(cache::NLNewtonCache, alg, i::Int)
  resize!(cache.ustep, i)
  resize!(cache.atmp, i)
  nothing
end

function nlsolve_cache_resize!(cache::NLNewtonConstantCache, alg, i::Int)
  nothing
end

function nlsolve_cache_resize!(cache::NLAndersonCache, alg, i::Int)
  resize!(cache.ustep, i)
  resize!(cache.atmp, i)
  resize!(cache.dzold, i)
  resize!(cache.z₊old, i)
  max_history = min(alg.nlsolve.max_history, alg.nlsolve.max_iter, i)
  prev_max_history = length(cache.Δz₊s)
  resize!(cache.γs, max_history)
  resize!(cache.Δz₊s, max_history)
  if max_history > prev_max_history
    for i in (max_history - prev_max_history):max_history
      cache.Δz₊s[i] = zero(z₊)
    end
  end
  cache.Q = typeof(cache.Q)(undef, i, max_history)
  cache.R = typeof(cache.R)(undef, max_history, max_history)
  nothing
end

function nlsolve_cache_resize!(cache::NLAndersonConstantCache, alg, i::Int)
  max_history = min(alg.nlsolve.max_history, alg.nlsolve.max_iter, i)
  resize!(cache.Δz₊s, max_history)
  cache.Q = typeof(cache.Q)(undef, i, max_history)
  cache.R = typeof(cache.R)(undef, max_history, max_history)
  resize!(cache.γs, max_history)
  nothing
end

function nlsolve_cache_resize!(cache::NLFunctionalCache, alg, i::Int)
  resize!(cache.ustep, i)
  resize!(cache.atmp, i)
  nothing
end

function nlsolve_cache_resize!(cache::NLFunctionalConstantCache, alg, i::Int)
  nothing
end
