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
  integrator.u_modified = false
end

function u_modified!(integrator::ODEIntegrator,bool::Bool)
  integrator.u_modified = bool
end

get_proposed_dt(integrator::ODEIntegrator) = integrator.dtpropose
set_proposed_dt!(integrator::ODEIntegrator,dt::Number) = (integrator.dtpropose = dt)

function set_proposed_dt!(integrator::ODEIntegrator,integrator2::ODEIntegrator)
  integrator.dtpropose = integrator2.dtpropose
  integrator.qold = integrator2.qold
  integrator.erracc = integrator2.erracc
  integrator.dtacc = integrator2.dtacc
end

@inline function DiffEqBase.get_du(integrator::ODEIntegrator)
  integrator.fsallast
end

@inline function DiffEqBase.get_du!(out,integrator::ODEIntegrator)
  out .= integrator.fsallast
end

#TODO: Bigger caches for most algorithms
@inline DiffEqBase.get_tmp_cache(integrator::ODEIntegrator) =
          get_tmp_cache(integrator::ODEIntegrator,integrator.alg,integrator.cache)
# avoid method ambiguity
for typ in (OrdinaryDiffEqAlgorithm,OrdinaryDiffEqNewtonAdaptiveAlgorithm,OrdinaryDiffEqRosenbrockAdaptiveAlgorithm)
  @eval @inline DiffEqBase.get_tmp_cache(integrator,alg::$typ,cache::OrdinaryDiffEqConstantCache) = nothing
end
@inline DiffEqBase.get_tmp_cache(integrator,alg::OrdinaryDiffEqAlgorithm,cache) = (cache.tmp,)
@inline DiffEqBase.get_tmp_cache(integrator,alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm,cache) = (cache.tmp,cache.atmp)
@inline DiffEqBase.get_tmp_cache(integrator,alg::OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,cache) = (cache.tmp,cache.linsolve_tmp)
@inline DiffEqBase.get_tmp_cache(integrator,alg::CompositeAlgorithm, cache) = get_tmp_cache(integrator, integrator.alg.algs[1], cache.caches[1])

cache_iter(integrator::ODEIntegrator) = _cache_iter(integrator.cache)

function add_tstop!(integrator::ODEIntegrator,t)
  integrator.tdir * (t - integrator.t) < 0 && error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
  push!(integrator.opts.tstops,t)
end

function DiffEqBase.add_saveat!(integrator::ODEIntegrator,t)
  integrator.tdir * (t - integrator.t) < 0 && error("Tried to add a saveat that is behind the current time. This is strictly forbidden")
  push!(integrator.opts.saveat,t)
end

resize!(integrator::ODEIntegrator,i::Int) = resize!(integrator,integrator.cache,i)
function resize!(integrator::ODEIntegrator,cache,i)
  for c in cache_iter(integrator)
    resize!(c,i)
  end
  resize_non_user_cache!(integrator,cache,i)
end

resize_non_user_cache!(integrator::ODEIntegrator,i::Int) = resize_non_user_cache!(integrator,integrator.cache,i)
deleteat_non_user_cache!(integrator::ODEIntegrator,i) = deleteat_non_user_cache!(integrator,integrator.cache,i)
addat_non_user_cache!(integrator::ODEIntegrator,i) = addat_non_user_cache!(integrator,integrator.cache,i)

resize_non_user_cache!(integrator::ODEIntegrator,cache,i) = nothing
function resize_non_user_cache!(integrator::ODEIntegrator,cache::RosenbrockMutableCache,i)
  cache.J = similar(cache.J,i,i)
  cache.W = similar(cache.W,i,i)
end

function resize_non_user_cache!(integrator::ODEIntegrator,cache::Union{GenericImplicitEulerCache,GenericTrapezoidCache},i)
  cache.nl_rhs = integrator.alg.nlsolve(Val{:init},cache.rhs,cache.u)
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
  for c in cache_iter(integrator)
    deleteat!(c,idxs)
  end
  deleteat_non_user_cache!(integrator,integrator.cache,idxs)
end

function addat!(integrator::ODEIntegrator,idxs)
  for c in cache_iter(integrator)
    addat!(c,idxs)
  end
  addat_non_user_cache!(integrator,integrator.cache,idxs)
end

function terminate!(integrator::ODEIntegrator)
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

  tstops_internal, saveat_internal, d_discontinuities_internal =
    tstop_saveat_disc_handling(tstops,saveat,d_discontinuities,
    integrator.tdir,(t0,tf),typeof(integrator.t))

  integrator.opts.tstops = tstops_internal
  integrator.opts.saveat = saveat_internal
  integrator.opts.d_discontinuities = d_discontinuities_internal

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
      if integrator.opts.save_idxs == nothing
        copyat_or_push!(integrator.sol.u,1,u0)
      else
        u_initial = u0[integrator.opts.save_idxs]
        copyat_or_push!(integrator.sol.u,1,u_initial,Val{false})
      end
    end
    if integrator.sol.u_analytic != nothing
      resize!(integrator.sol.u_analytic,0)
    end
    if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
      resize!(integrator.sol.alg_choice,resize_start)
    end
    integrator.saveiter = resize_start
  end
  integrator.iter = 0
  integrator.success_iter = 0
  integrator.u_modified = false

  # full re-initialize the PI in timestepping
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
