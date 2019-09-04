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

build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c,iip) =
  build_nlsolver(alg,alg.nlsolve,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c,iip)

function build_nlsolver(alg,nlalg::Union{NLFunctional,NLAnderson,NLNewton},u,uprev,p,t,dt,
                        f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c,::Val{true})
  @unpack κ, fast_convergence_cutoff = nlalg

  # define additional fields of cache of non-linear solver
  z = similar(u); dz = similar(u); tmp = similar(u); b = similar(u)
  k = zero(rate_prototype); ustep = similar(u)
  tstep = zero(t)
  atmp = similar(u, uEltypeNoUnits)

  uTolType = real(uBottomEltypeNoUnits)

  # create cache of non-linear solver
  if nlalg isa NLNewton
    nf = nlsolve_f(f, alg)

    if islinear(f)
      du1 = rate_prototype
      uf = nothing
      jac_config = nothing
      linsolve = alg.linsolve(Val{:init},nf,u)
    else
      du1 = zero(rate_prototype)
      uf = build_uf(alg,nf,t,p,Val(true))
      jac_config = build_jac_config(alg,nf,uf,du1,uprev,u,tmp,dz)
      linsolve = alg.linsolve(Val{:init},uf,u)
    end

    # TODO: check if the solver is iterative
    weight = similar(u)

    tType = typeof(t)
    invγdt = inv(oneunit(t) * one(uTolType))

    J, W = build_J_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits,Val(true))

    nlcache = NLNewtonCache(ustep,tstep,atmp,dz,J,W,true,tType(dt),du1,uf,jac_config,
                            linsolve,weight,invγdt,tType(nlalg.new_W_dt_cutoff))
  elseif nlalg isa NLFunctional
    nlcache = NLFunctionalCache(ustep,tstep,atmp,dz)
  elseif nlalg isa NLAnderson
    max_history = min(nlalg.max_history, nlalg.max_iter, length(z))
    Δz₊s = [zero(z) for i in 1:max_history]
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)

    dzold = zero(z)
    z₊old = zero(z)

    nlcache = NLAndersonCache(ustep,tstep,atmp,dz,dzold,z₊old,Δz₊s,Q,R,γs,0,nlalg.aa_start,
                              nlalg.droptol)
  end

  # create non-linear solver
  NLSolver{typeof(nlalg),true,typeof(z),typeof(k),uTolType,typeof(κ),typeof(γ),typeof(c),
           typeof(fast_convergence_cutoff),typeof(nlcache)}(
             z,tmp,b,k,nlalg,one(uTolType),κ,γ,c,nlalg.max_iter,10000,Convergence,
             fast_convergence_cutoff,nlcache)
end

function build_nlsolver(alg,nlalg::Union{NLFunctional,NLAnderson,NLNewton},u,uprev,p,t,dt,
                        f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c,::Val{false})
  @unpack κ, fast_convergence_cutoff = nlalg

  # define additional fields of cache of non-linear solver (all aliased)
  z = uprev; tmp = z; b = z; k = rate_prototype
  tstep = zero(t)

  uTolType = real(uBottomEltypeNoUnits)

  # create cache of non-linear solver
  if nlalg isa NLNewton
    nf = nlsolve_f(f, alg)
    uf = build_uf(alg,nf,t,p,Val(false))

    invγdt = inv(oneunit(t) * one(uTolType))

    J, W = build_J_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits,Val(false))

    nlcache = NLNewtonConstantCache(tstep,J,W,uf,invγdt,typeof(t)(nlalg.new_W_dt_cutoff))
  elseif nlalg isa NLFunctional
    nlcache = NLFunctionalConstantCache(tstep)
  elseif nlalg isa NLAnderson
    max_history = min(nlalg.max_history, nlalg.max_iter, length(z))
    Δz₊s = Vector{typeof(z)}(undef, max_history)
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)

    dz = u
    dzold = u
    z₊old = u

    nlcache = NLAndersonConstantCache(tstep,dz,dzold,z₊old,Δz₊s,Q,R,γs,0,nlalg.aa_start,nlalg.droptol)
  end

  # create non-linear solver
  NLSolver{typeof(nlalg),false,typeof(z),typeof(k),uTolType,typeof(κ),typeof(γ),typeof(c),
           typeof(fast_convergence_cutoff),typeof(nlcache)}(
             z,tmp,b,k,nlalg,one(uTolType),κ,γ,c,nlalg.max_iter,10000,Convergence,
             fast_convergence_cutoff,nlcache)
end

## Anderson acceleration

"""
    anderson(z, cache, integrator)

Return the next iterate of the fixed-point iteration `z = g(z)` by performing Anderson
acceleration based on the current iterate `z` and the settings and history in the `cache`.
"""
@muladd function anderson(z, cache, integrator)
  @unpack dz,Δz₊s,z₊old,dzold,R,Q,γs,history,droptol = cache

  # increase size of history
  history += 1

  # remove oldest history if maximum size is exceeded
  max_history = length(Δz₊s)
  if history > max_history
    # circularly shift differences of z₊
    for i in 1:(max_history-1)
      Δz₊s[i] = Δz₊s[i + 1]
    end

    # delete left-most column of QR decomposition
    qrdelete!(Q, R, max_history)

    # update size of history
    history = max_history
  end

  # update history of differences of z₊
  Δz₊s[history] = @.. z - z₊old

  # replace/add difference of residuals as right-most column to QR decomposition
  qradd!(Q, R, _vec(dz .- dzold), history)

  # update cached values
  cache.dzold = dz
  cache.z₊old = z

  # define current Q and R matrices
  Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

  # check condition (TODO: incremental estimation)
  if droptol !== nothing
    while cond(R) > droptol && history > 1
      qrdelete!(Q, R, history)
      history -= 1
      Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
    end
  end

  # save updated history
  cache.history = history

  # solve least squares problem
  γscur = view(γs, 1:history)
  ldiv!(Rcur, mul!(γscur, Qcur', _vec(dz)))
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nsolve += 1
  end

  # update next iterate
  for i in 1:history
    z = @.. z - γs[i] * Δz₊s[i]
  end

  z
end

"""
    anderson!(z, cache, integrator)

Update the current iterate `z` of the fixed-point iteration `z = g(z)` in-place
by performing Anderson acceleration based on the settings and history in the `cache`.
"""
@muladd function anderson!(z, cache, integrator)
  @unpack dz,z₊old,dzold,Δz₊s,γs,R,Q,history,droptol = cache

  # increase size of history
  history += 1

  # remove oldest history if maximum size is exceeded
  max_history = length(Δz₊s)
  if history > max_history
    # circularly shift differences of z₊
    ptr = Δz₊s[1]
    for i in 1:(max_history-1)
      Δz₊s[i] = Δz₊s[i + 1]
    end
    Δz₊s[max_history] = ptr

    # delete left-most column of QR decomposition
    qrdelete!(Q, R, max_history)

    # update size of history
    history = max_history
  end

  # update history of differences of z₊
  @.. Δz₊s[history] = z - z₊old

  # replace/add difference of residuals as right-most column to QR decomposition
  @.. dzold = dz - dzold
  qradd!(Q, R, vec(dzold), history)

  # update cached values
  @.. dzold = dz
  @.. z₊old = z

  # define current Q and R matrices
  Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

  # check condition (TODO: incremental estimation)
  if droptol !== nothing
    while cond(R) > droptol && history > 1
      qrdelete!(Q, R, history)
      history -= 1
      Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
    end
  end

  # save updated history
  cache.history = history

  # solve least squares problem
  γscur = view(γs, 1:history)
  ldiv!(Rcur, mul!(γscur, Qcur', vec(dz)))
  if DiffEqBase.has_destats(integrator)
    integrator.destats.nsolve += 1
  end

  # update next iterate
  for i in 1:history
    @.. z = z - γs[i] * Δz₊s[i]
  end

  nothing
end  

function nlsolve_resize!(integrator::DiffEqBase.DEIntegrator, i::Int)
  if !isdefined(integrator.cache, :nlsolver)
    return nothing
  end
  alg = integrator.alg; nlsolver = integrator.cache.nlsolver
  if nlsolver isa AbstractArray
    for idx in eachindex(nlsolver) # looping because we may have multiple nlsolver for threaded case
      _nlsolver = nlsolver[idx]
      @unpack z,tmp,ztmp,k,cache = _nlsolver
      # doubt: if these fields are always going to be in alg cache too, then we shouldnt do this here.
      # double resize doesn't do any bad I think though
      resize!(z,i)
      resize!(tmp,i)
      resize!(ztmp,i)
      resize!(k,i)
      nlsolve_cache_resize!(cache,alg,i)
    end
  else
    @unpack z,tmp,ztmp,k,cache = nlsolver
    resize!(z,i)
    resize!(tmp,i)
    resize!(ztmp,i)
    resize!(k,i)
    nlsolve_cache_resize!(cache,alg,i)
  end
  nothing
end

function nlsolve_cache_resize!(cache::NLNewtonCache, alg, i::Int)
  resize!(cache.ustep, i)
  resize!(cache.atmp, i)
  resize!(cache.dz,i)
  resize!(cache.du1, i)
  if cache.jac_config !== nothing
    resize_jac_config!(cache.jac_config, i)
  end
  resize!(cache.weight, i)

  nothing
end

function nlsolve_cache_resize!(cache::NLNewtonConstantCache, alg, i::Int)
  nothing
end

function nlsolve_cache_resize!(cache::NLAndersonCache, alg, i::Int)
  resize!(cache.ustep, i)
  resize!(cache.atmp, i)
  resize!(cache.dz, i)
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
  resize!(cache.dz, i)
  nothing
end

function nlsolve_cache_resize!(cache::NLFunctionalConstantCache, alg, i::Int)
  nothing
end
