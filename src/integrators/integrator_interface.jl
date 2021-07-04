# We want to make sure that the first argument of change_t_via_interpolation!
# is specialized, yet, it needs to take both ODEIntegrator and DDEIntegrator.
# Hence, we need to have two separate functions.
function _change_t_via_interpolation!(integrator,t,modify_save_endpoint::Type{Val{T}}) where T
  # Can get rid of an allocation here with a function
  # get_tmp_arr(integrator.cache) which gives a pointer to some
  # cache array which can be modified.
  if integrator.tdir*t < integrator.tdir*integrator.tprev
    error("Current interpolant only works between tprev and t")
  elseif t != integrator.t
    if is_constant_cache(integrator.cache)
      integrator.u = integrator(t)
    else
      integrator(integrator.u,t)
    end
    integrator.t = t
    integrator.dt = integrator.t - integrator.tprev
    DiffEqBase.reeval_internals_due_to_modification!(integrator)
    if T
      solution_endpoint_match_cur_integrator!(integrator)
    end
  end
end
DiffEqBase.change_t_via_interpolation!(integrator::ODEIntegrator,
                                        t,modify_save_endpoint::Type{Val{T}}=Val{false}) where T =
                                          _change_t_via_interpolation!(integrator,t,modify_save_endpoint)

function DiffEqBase.reeval_internals_due_to_modification!(integrator::ODEIntegrator)
  if integrator.opts.calck
    resize!(integrator.k,integrator.kshortsize) # Reset k for next step!
    alg = unwrap_alg(integrator, false)
    if typeof(alg) <: BS5 || typeof(alg) <: Vern6 || typeof(alg) <: Vern7 ||
       typeof(alg) <: Vern8 || typeof(alg) <: Vern9
       DiffEqBase.addsteps!(integrator,integrator.f,true,false,!alg.lazy)
    else
      DiffEqBase.addsteps!(integrator,integrator.f,true,false)
    end
  end

  if integrator.isdae
    DiffEqBase.initialize_dae!(integrator)
  end

  integrator.u_modified = false
end

function u_modified!(integrator::ODEIntegrator,bool::Bool)
  integrator.u_modified = bool
end

get_proposed_dt(integrator::ODEIntegrator) = ifelse(integrator.opts.adaptive, integrator.dtpropose, integrator.dtcache)
set_proposed_dt!(integrator::ODEIntegrator,dt::Number) = (integrator.dtpropose = dt; integrator.dtcache = dt)

function set_proposed_dt!(integrator::ODEIntegrator,integrator2::ODEIntegrator)
  integrator.dtpropose = integrator2.dtpropose
  integrator.dtcache = integrator2.dtcache
  integrator.qold = integrator2.qold
  integrator.erracc = integrator2.erracc
  integrator.dtacc = integrator2.dtacc
end

@inline function DiffEqBase.get_du(integrator::ODEIntegrator)
  return if isdefined(integrator, :fsallast)
    integrator.fsallast
  else
    integrator(integrator.t, Val{1})
  end
end

@inline function DiffEqBase.get_du!(out,integrator::ODEIntegrator)
  if typeof(integrator.cache) <: FunctionMapCache
    out .= integrator.cache.tmp
  else
    return if isdefined(integrator, :fsallast)
      out .= integrator.fsallast
    else
      integrator(out, integrator.t, Val{1})
    end
  end
end

#TODO: Bigger caches for most algorithms
@inline DiffEqBase.get_tmp_cache(integrator::ODEIntegrator) =
          get_tmp_cache(integrator::ODEIntegrator,integrator.alg,integrator.cache)
# avoid method ambiguity
for typ in (OrdinaryDiffEqAlgorithm,Union{RadauIIA3,RadauIIA5},OrdinaryDiffEqNewtonAdaptiveAlgorithm,OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,Union{SSPRK22,SSPRK33,SSPRK53_2N1,SSPRK53_2N2,SSPRK43,SSPRK432,SSPRK932})
  @eval @inline DiffEqBase.get_tmp_cache(integrator,alg::$typ,cache::OrdinaryDiffEqConstantCache) = nothing
end

# the ordering of the cache arrays is important!!!
@inline DiffEqBase.get_tmp_cache(integrator,alg::OrdinaryDiffEqAlgorithm,cache) = (cache.tmp,)
@inline DiffEqBase.get_tmp_cache(integrator,alg::Union{RadauIIA3,RadauIIA5},cache) = (cache.tmp,cache.atmp)
@inline DiffEqBase.get_tmp_cache(integrator,alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm,cache) = (cache.nlsolver.tmp,cache.atmp)
@inline DiffEqBase.get_tmp_cache(integrator,alg::OrdinaryDiffEqNewtonAlgorithm,cache) = (cache.nlsolver.tmp,cache.nlsolver.z)
@inline DiffEqBase.get_tmp_cache(integrator,alg::OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,cache) = (cache.tmp,cache.linsolve_tmp)
@inline DiffEqBase.get_tmp_cache(integrator,alg::Union{SSPRK22,SSPRK33,SSPRK53_2N1,SSPRK53_2N2,SSPRK43,SSPRK432,SSPRK932},cache) = (cache.k,)
@inline DiffEqBase.get_tmp_cache(integrator,alg::OrdinaryDiffEqImplicitExtrapolationAlgorithm,cache) = (cache.tmp,cache.utilde)
@inline DiffEqBase.get_tmp_cache(integrator,alg::OrdinaryDiffEqAdaptiveExponentialAlgorithm,cache) = (cache.tmp,cache.utilde)
@inline DiffEqBase.get_tmp_cache(integrator,alg::OrdinaryDiffEqExponentialAlgorithm,cache) = (cache.tmp,cache.dz)
@inline DiffEqBase.get_tmp_cache(integrator,alg::OrdinaryDiffEqLinearExponentialAlgorithm,cache) = (cache.tmp,)
@inline DiffEqBase.get_tmp_cache(integrator,alg::CompositeAlgorithm, cache) = get_tmp_cache(integrator, integrator.alg.algs[1], cache.caches[1])
@inline DiffEqBase.get_tmp_cache(integrator,alg::DAEAlgorithm,cache) = (cache.nlsolver.cache.dz,cache.atmp)

full_cache(integrator::ODEIntegrator) = full_cache(integrator.cache)
full_cache(integrator::CompositeCache) = Iterators.flatten(full_cache(c) for c in integrator.caches)

function add_tstop!(integrator::ODEIntegrator,t)
  integrator.tdir * (t - integrator.t) < zero(integrator.t) && error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
  push!(integrator.opts.tstops, integrator.tdir * t)
end

function DiffEqBase.add_saveat!(integrator::ODEIntegrator,t)
  integrator.tdir * (t - integrator.t) < zero(integrator.t) && error("Tried to add a saveat that is behind the current time. This is strictly forbidden")
  push!(integrator.opts.saveat, integrator.tdir * t)
end

function resize!(integrator::ODEIntegrator, i::Int)
  @unpack cache = integrator

  for c in full_cache(cache)
    # Skip nothings which may exist in the cache since extra variables
    # may be required for things like units
    c !== nothing && resize!(c,i)
  end
  resize_nlsolver!(integrator, i)
  resize_J_W!(cache, integrator, i)
  resize_non_user_cache!(integrator, cache, i)
end
# we can't use resize!(..., i::Union{Int, NTuple{N,Int}}) where {N} because of method ambiguities with DiffEqBase
function resize!(integrator::ODEIntegrator, i::NTuple{N,Int}) where {N}
  @unpack cache = integrator

  for c in full_cache(cache)
    resize!(c,i)
  end
  # TODO the parts below need to be adapted for implicit methods
  isdefined(integrator.cache, :nlsolver) && resize_nlsolver!(integrator, i)
  resize_J_W!(cache, integrator, i)
  resize_non_user_cache!(integrator, cache, i)
end

function resize_J_W!(cache, integrator, i)
  (isdefined(cache, :J) && isdefined(cache, :W)) || return

  @unpack f = integrator

  if cache.W isa WOperator
    nf = nlsolve_f(f, integrator.alg)
    islin = f isa Union{ODEFunction,SplitFunction} && islinear(nf.f)
    if !islin
      if isa(cache.J, DiffEqBase.AbstractDiffEqLinearOperator)
        resize!(cache.J,i)
      else
        J = similar(f.jac_prototype, i, i)
        J = DiffEqArrayOperator(J; update_func=f.jac)
      end
      cache.W = WOperator{DiffEqBase.isinplace(integrator.sol.prob)}(
                          f.mass_matrix, integrator.dt, cache.J, integrator.u;
                          transform = cache.W.transform)
      cache.J = cache.W.J
    end
  else
    if cache.J !== nothing
      cache.J = similar(cache.J, i, i)
    end
    cache.W = similar(cache.W, i, i)
  end

  nothing
end

resize_non_user_cache!(integrator::ODEIntegrator,i::Int) = resize_non_user_cache!(integrator,integrator.cache,i)
deleteat_non_user_cache!(integrator::ODEIntegrator,i) = deleteat_non_user_cache!(integrator,integrator.cache,i)
addat_non_user_cache!(integrator::ODEIntegrator,i) = addat_non_user_cache!(integrator,integrator.cache,i)

resize_non_user_cache!(integrator::ODEIntegrator,cache,i) = nothing

function resize_non_user_cache!(integrator::ODEIntegrator,cache::CompositeCache,i)
  for _cache in cache.caches
    resize_non_user_cache!(integrator,_cache,i)
  end
end

function deleteat_non_user_cache!(integrator::ODEIntegrator,cache::CompositeCache,i)
  for _cache in cache.caches
    deleteat_non_user_cache!(integrator,_cache,i)
  end
end

function addat_non_user_cache!(integrator::ODEIntegrator,cache::CompositeCache,i)
  for _cache in cache.caches
    addat_non_user_cache!(integrator,_cache,i)
  end
end

function resize_non_user_cache!(integrator::ODEIntegrator,
                      cache::RosenbrockMutableCache,i)
  cache.J = similar(cache.J,i,i)
  cache.W = similar(cache.W,i,i)
  resize_jac_config!(cache.jac_config, i)
  resize_grad_config!(cache.grad_config, i)
  nothing
end

function deleteat_non_user_cache!(integrator::ODEIntegrator,cache,idxs)
  # ordering doesn't matter in deterministic cache, so just resize
  # to match the size of u
  i = length(integrator.u)
  resize_non_user_cache!(integrator,cache,i)
end

function addat_non_user_cache!(integrator::ODEIntegrator,cache,idxs)
  # ordering doesn't matter in deterministic cache, so just resize
  # to match the size of u
  i = length(integrator.u)
  resize_non_user_cache!(integrator,cache,i)
end

function deleteat!(integrator::ODEIntegrator,idxs)
  for c in full_cache(integrator)
    deleteat!(c,idxs)
  end
  deleteat_non_user_cache!(integrator,integrator.cache,idxs)
end

function addat!(integrator::ODEIntegrator,idxs)
  for c in full_cache(integrator)
    addat!(c,idxs)
  end
  addat_non_user_cache!(integrator,integrator.cache,idxs)
end

function terminate!(integrator::ODEIntegrator, retcode = :Terminated)
  integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, retcode)
  integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
end

DiffEqBase.has_reinit(integrator::ODEIntegrator) = true
function DiffEqBase.reinit!(integrator::ODEIntegrator,u0 = integrator.sol.prob.u0;
  t0 = integrator.sol.prob.tspan[1], tf = integrator.sol.prob.tspan[2],
  erase_sol = true,
  tstops = integrator.opts.tstops_cache,
  saveat = integrator.opts.saveat_cache,
  d_discontinuities = integrator.opts.d_discontinuities_cache,
  reset_dt = (integrator.dtcache == zero(integrator.dt)) && integrator.opts.adaptive,
  reinit_callbacks = true, initialize_save = true,
  reinit_cache = true)

  if isinplace(integrator.sol.prob)
    recursivecopy!(integrator.u,u0)
    recursivecopy!(integrator.uprev,integrator.u)
  else
    integrator.u = u0
    integrator.uprev = integrator.u
  end

  if alg_extrapolates(integrator.alg)
    if isinplace(integrator.sol.prob)
      recursivecopy!(integrator.uprev2,integrator.uprev)
    else
      integrator.uprev2 = integrator.uprev
    end
  end

  integrator.t = t0
  integrator.tprev = t0

  tType = typeof(integrator.t)
  tspan = (tType(t0), tType(tf))
  integrator.opts.tstops = initialize_tstops(tType, tstops, d_discontinuities, tspan)
  integrator.opts.saveat = initialize_saveat(tType, saveat, tspan)
  integrator.opts.d_discontinuities =
    initialize_d_discontinuities(tType, d_discontinuities, tspan)

  if erase_sol
    if integrator.opts.save_start
      resize_start = 1
    else
      resize_start = 0
    end
    resize!(integrator.sol.u,resize_start)
    resize!(integrator.sol.t,resize_start)
    resize!(integrator.sol.k,resize_start)
    if integrator.opts.save_start
      copyat_or_push!(integrator.sol.t,1,t0)
      if integrator.opts.save_idxs === nothing
        copyat_or_push!(integrator.sol.u,1,u0)
      else
        u_initial = u0[integrator.opts.save_idxs]
        copyat_or_push!(integrator.sol.u,1,u_initial,Val{false})
      end
    end
    if integrator.sol.u_analytic !== nothing
      resize!(integrator.sol.u_analytic,0)
    end
    if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
      resize!(integrator.sol.alg_choice,resize_start)
    end
    integrator.saveiter = resize_start
    if integrator.opts.dense
      integrator.saveiter_dense = resize_start
    end
  end
  integrator.iter = 0
  integrator.success_iter = 0
  integrator.u_modified = false

  # full re-initialize the PI in timestepping
  reinit!(integrator, integrator.opts.controller)
  integrator.qold = integrator.opts.qoldinit
  integrator.q11 = typeof(integrator.q11)(1)
  integrator.erracc = typeof(integrator.erracc)(1)
  integrator.dtacc = typeof(integrator.dtacc)(1)

  if reset_dt
    auto_dt_reset!(integrator)
  end

  if reinit_callbacks
    initialize_callbacks!(integrator, initialize_save)
  end

  if reinit_cache
    initialize!(integrator,integrator.cache)
  end
end

function DiffEqBase.auto_dt_reset!(integrator::ODEIntegrator)
  integrator.dt = ode_determine_initdt(integrator.u,integrator.t,
  integrator.tdir,integrator.opts.dtmax,integrator.opts.abstol,integrator.opts.reltol,
  integrator.opts.internalnorm,integrator.sol.prob,integrator)
  integrator.destats.nf += 2
end

function DiffEqBase.set_t!(integrator::ODEIntegrator, t::Real)
  if integrator.opts.save_everystep
    error("Integrator time cannot be reset unless it is initialized",
          " with save_everystep=false")
  end
  if alg_extrapolates(integrator.alg) || !isdtchangeable(integrator.alg)
    reinit!(integrator, integrator.u;
            t0 = t,
            reset_dt = false,
            reinit_callbacks = false,
            reinit_cache = false)
  else
    integrator.t = t
  end
end

function DiffEqBase.set_u!(integrator::ODEIntegrator, u)
  if integrator.opts.save_everystep
    error("Integrator state cannot be reset unless it is initialized",
          " with save_everystep=false")
  end
  integrator.u = u
  u_modified!(integrator, true)
end

DiffEqBase.has_destats(i::ODEIntegrator) = true
