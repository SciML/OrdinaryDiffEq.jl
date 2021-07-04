get_status(nlsolver::AbstractNLSolver) = nlsolver.status
get_new_W_γdt_cutoff(nlsolver::AbstractNLSolver) = nlsolver.cache.new_W_γdt_cutoff
# handle FIRK
get_new_W_γdt_cutoff(alg::NewtonAlgorithm) = alg.new_W_γdt_cutoff

nlsolvefail(nlsolver::AbstractNLSolver) = nlsolvefail(get_status(nlsolver))
nlsolvefail(status::NLStatus) = Int8(status) <= 0

isnewton(::Any) = false
isnewton(nlsolver::AbstractNLSolver) = isnewton(nlsolver.cache)
isnewton(::AbstractNLSolverCache) = false
isnewton(::Union{NLNewtonCache,NLNewtonConstantCache}) = true

isJcurrent(nlsolver::AbstractNLSolver, integrator) = integrator.t == nlsolver.cache.J_t
isfirstcall(nlsolver::AbstractNLSolver) = nlsolver.cache.firstcall
isfirststage(nlsolver::AbstractNLSolver) = nlsolver.cache.firststage
setfirststage!(nlsolver::AbstractNLSolver, val::Bool) = setfirststage!(nlsolver.cache, val)
setfirststage!(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}, val::Bool) = (nlcache.firststage = val)
setfirststage!(::Any, val::Bool) = nothing
markfirststage!(nlsolver::AbstractNLSolver) = setfirststage!(nlsolver, true)

set_new_W!(nlsolver::AbstractNLSolver, val::Bool)::Bool = set_new_W!(nlsolver.cache, val)
set_new_W!(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}, val::Bool)::Bool =
  nlcache.new_W = val
get_new_W!(nlsolver::AbstractNLSolver)::Bool = get_new_W!(nlsolver.cache)
get_new_W!(nlcache::Union{NLNewtonCache,NLNewtonConstantCache})::Bool = nlcache.new_W
get_new_W!(::AbstractNLSolverCache)::Bool = true

get_W(nlsolver::AbstractNLSolver) = get_W(nlsolver.cache)
get_W(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}) = nlcache.W

set_W_γdt!(nlsolver::AbstractNLSolver, W_γdt) = set_W_γdt!(nlsolver.cache, W_γdt)
function set_W_γdt!(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}, W_γdt)
  nlcache.W_γdt = W_γdt
  W_γdt
end

du_cache(nlsolver::AbstractNLSolver) = du_cache(nlsolver.cache)
du_cache(::AbstractNLSolverCache) = nothing
du_cache(nlcache::Union{NLFunctionalCache,NLAndersonCache,NLNewtonCache}) = (nlcache.k,)

function du_alias_or_new(nlsolver::AbstractNLSolver, rate_prototype)
  _du_cache = du_cache(nlsolver)
  if _du_cache === nothing
    zero(rate_prototype)
  else
    first(_du_cache)
  end
end

mutable struct DAEResidualJacobianWrapper{AD,F,pType,duType,uType,alphaType,gammaType,tmpType,uprevType,tType} <: Function
  f::F
  p::pType
  tmp_du::duType
  tmp_u::uType
  α::alphaType
  invγdt::gammaType
  tmp::tmpType
  uprev::uprevType
  t::tType
  function DAEResidualJacobianWrapper(alg,f,p,α,invγdt,tmp,uprev,t)
    isautodiff = alg_autodiff(alg)
    if isautodiff
      tmp_du = DiffEqBase.dualcache(uprev)
      tmp_u = DiffEqBase.dualcache(uprev)
    else
      tmp_du = similar(uprev)
      tmp_u = similar(uprev)
    end
    new{isautodiff,typeof(f),typeof(p),typeof(tmp_du),typeof(tmp_u),typeof(α),typeof(invγdt),typeof(tmp),typeof(uprev),typeof(t)}(f,p,tmp_du,tmp_u,α,invγdt,tmp,uprev,t)
  end
end

is_autodiff(m::DAEResidualJacobianWrapper{AD}) where AD = AD

function (m::DAEResidualJacobianWrapper)(out,x)
  if is_autodiff(m)
    tmp_du = DiffEqBase.get_tmp(m.tmp_du, x)
    tmp_u = DiffEqBase.get_tmp(m.tmp_u, x)
  else
    tmp_du = m.tmp_du
    tmp_u = m.tmp_u
  end
  @. tmp_du = (m.α * x + m.tmp) * m.invγdt
  @. tmp_u = x + m.uprev
  m.f(out, tmp_du, tmp_u, m.p, m.t)
end

mutable struct DAEResidualDerivativeWrapper{F,pType,alphaType,gammaType,tmpType,uprevType,tType} <: Function
  f::F
  p::pType
  α::alphaType
  invγdt::gammaType
  tmp::tmpType
  uprev::uprevType
  t::tType
end

function (m::DAEResidualDerivativeWrapper)(x)
  tmp_du = (m.α * x + m.tmp) * m.invγdt
  tmp_u = x + m.uprev
  m.f(tmp_du, tmp_u, m.p, m.t)
end

DiffEqBase.has_jac(f::DAEResidualJacobianWrapper) = DiffEqBase.has_jac(f.f)
DiffEqBase.has_Wfact(f::DAEResidualJacobianWrapper) = DiffEqBase.has_Wfact(f.f)
DiffEqBase.has_Wfact_t(f::DAEResidualJacobianWrapper) = DiffEqBase.has_Wfact_t(f.f)

DiffEqBase.has_jac(f::DAEResidualDerivativeWrapper) = DiffEqBase.has_jac(f.f)
DiffEqBase.has_Wfact(f::DAEResidualDerivativeWrapper) = DiffEqBase.has_Wfact(f.f)
DiffEqBase.has_Wfact_t(f::DAEResidualDerivativeWrapper) = DiffEqBase.has_Wfact_t(f.f)

function build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                        tTypeNoUnits,γ,c,iip)
  build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                        tTypeNoUnits,γ,c,1,iip)
end

function build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                        tTypeNoUnits,γ,c,α,iip)
  build_nlsolver(alg,alg.nlsolve,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,
                 uBottomEltypeNoUnits,tTypeNoUnits,γ,c,α,iip)
end

function build_nlsolver(alg,nlalg::Union{NLFunctional,NLAnderson,NLNewton},u,uprev,p,t,dt,
                        f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                        γ,c,α,::Val{true})
  #TODO
  #nlalg = DiffEqBase.handle_defaults(alg, nlalg)
  # define unitless type
  uTolType = real(uBottomEltypeNoUnits)
  isdae = alg isa DAEAlgorithm

  # define fields of non-linear solver
  z = zero(u); tmp = zero(u); ztmp = zero(u)

  # build cache of non-linear solver
  ustep = zero(u)
  tstep = zero(t)
  k = zero(rate_prototype)
  atmp = similar(u, uEltypeNoUnits)
  dz = zero(u)

  if nlalg isa NLNewton
    nf = nlsolve_f(f, alg)

    if islinear(f)
      du1 = rate_prototype
      uf = nothing
      jac_config = nothing
      linsolve = alg.linsolve(Val{:init},nf,u)
    else
      du1 = zero(rate_prototype)
      if isdae
        uf = DAEResidualJacobianWrapper(alg,f,p,α,inv(γ*dt),tmp,uprev,t)
      else
        uf = build_uf(alg,nf,t,p,Val(true))
      end
      jac_config = build_jac_config(alg,nf,uf,du1,uprev,u,ztmp,dz)
      linsolve = alg.linsolve(Val{:init},uf,u)
    end

    # TODO: check if the solver is iterative
    weight = zero(u)

    tType = typeof(t)
    invγdt = inv(oneunit(t) * one(uTolType))

    J, W = build_J_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits,Val(true))

    nlcache = NLNewtonCache(ustep,tstep,k,atmp,dz,J,W,true,true,true,tType(dt),du1,uf,jac_config,
                            linsolve,weight,invγdt,tType(nlalg.new_W_dt_cutoff),t)
  elseif nlalg isa NLFunctional
    nlcache = NLFunctionalCache(ustep,tstep,k,atmp,dz)
  elseif nlalg isa NLAnderson
    max_history = min(nlalg.max_history, nlalg.max_iter, length(z))
    Δz₊s = [zero(z) for i in 1:max_history]
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)

    dzold = zero(z)
    z₊old = zero(z)

    nlcache = NLAndersonCache(ustep,tstep,atmp,k,dz,dzold,z₊old,Δz₊s,Q,R,γs,0,
                              nlalg.aa_start,nlalg.droptol)
  end

  # build non-linear solver
  ηold = one(t)

  NLSolver{true,tTypeNoUnits}(
    z,tmp,ztmp,γ,c,α,nlalg,nlalg.κ,
    nlalg.fast_convergence_cutoff,ηold,0,nlalg.max_iter,Divergence,
    nlcache)
end

function build_nlsolver(alg,nlalg::Union{NLFunctional,NLAnderson,NLNewton},u,uprev,p,t,dt,
                        f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                        γ,c,α,::Val{false})
  #TODO
  #nlalg = DiffEqBase.handle_defaults(alg, nlalg)
  # define unitless type
  uTolType = real(uBottomEltypeNoUnits)
  isdae = alg isa DAEAlgorithm

  # define fields of non-linear solver
  z = u; tmp = u; ztmp = u

  # build cache of non-linear solver
  tstep = zero(t)

  if nlalg isa NLNewton
    nf = nlsolve_f(f, alg)
    if isdae
      uf = DAEResidualDerivativeWrapper(f,p,α,inv(γ*dt),tmp,uprev,t)
    else
      uf = build_uf(alg,nf,t,p,Val(false))
    end

    tType = typeof(t)
    invγdt = inv(oneunit(t) * one(uTolType))

    J, W = build_J_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits,Val(false))

    nlcache = NLNewtonConstantCache(tstep,J,W,true,true,true,tType(dt),uf,invγdt,tType(nlalg.new_W_dt_cutoff),t)
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

  # build non-linear solver
  ηold = one(tTypeNoUnits)
  NLSolver{false,tTypeNoUnits}(
    z, tmp, ztmp, γ, c, α, nlalg, nlalg.κ,
    nlalg.fast_convergence_cutoff, ηold, 0, nlalg.max_iter, Divergence,
    nlcache)
end

## Anderson acceleration

"""
    anderson(z, cache)

Return the next iterate of the fixed-point iteration `z = g(z)` by performing Anderson
acceleration based on the current iterate `z` and the settings and history in the `cache`.
"""
@muladd function anderson(z, cache)
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

  # update next iterate
  for i in 1:history
    z = @.. z - γs[i] * Δz₊s[i]
  end

  z
end

"""
    anderson!(z, cache)

Update the current iterate `z` of the fixed-point iteration `z = g(z)` in-place
by performing Anderson acceleration based on the settings and history in the `cache`.
"""
@muladd function anderson!(z, cache)
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

  # update next iterate
  for i in 1:history
    @.. z = z - γs[i] * Δz₊s[i]
  end

  nothing
end

## resize

function resize_nlsolver!(integrator::DiffEqBase.DEIntegrator, i::Int)
  isdefined(integrator.cache, :nlsolver) || return

  @unpack nlsolver = integrator.cache

  if nlsolver isa AbstractArray
    for idx in eachindex(nlsolver)
      resize!(nlsolver[idx], integrator, i)
    end
  else
    resize!(nlsolver, integrator, i)
  end

  nlsolver.alg isa NLNewton && resize!(nlsolver.cache.linsolve,i)

  # make it reset everything since the caches changed size!
  nlsolver.cache.firstcall = true

  nothing
end

function Base.resize!(nlsolver::AbstractNLSolver, integrator, i::Int)
  resize!(nlsolver.z, i)
  resize!(nlsolver.tmp, i)
  resize!(nlsolver.ztmp, i)

  resize!(nlsolver.cache, nlsolver, integrator, i)
end

## default: dispatch only on the cache
Base.resize!(cache::AbstractNLSolverCache, nlsolver, integrator, i::Int) =
  Base.resize!(cache, i)
