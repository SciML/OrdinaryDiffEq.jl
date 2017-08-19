function change_t_via_interpolation!{T}(integrator,t,modify_save_endpoint::Type{Val{T}}=Val{false})
  # Can get rid of an allocation here with a function
  # get_tmp_arr(integrator.cache) which gives a pointer to some
  # cache array which can be modified.
  if integrator.tdir*t < integrator.tdir*integrator.tprev
    error("Current interpolant only works between tprev and t")
  elseif t != integrator.t

    if typeof(integrator.u) <: AbstractArray
      integrator(integrator.u,t)
    else
      integrator.u = integrator(t)
    end
    integrator.t = t
    integrator.dt = integrator.t - integrator.tprev
    reeval_internals_due_to_modification!(integrator)
    if T
      solution_endpoint_match_cur_integrator!(integrator)
    end
  end
end

function reeval_internals_due_to_modification!(integrator)
  if integrator.opts.calck
    resize!(integrator.k,integrator.kshortsize) # Reset k for next step!
    ode_addsteps!(integrator,integrator.f,Val{true},Val{false})
  end
  integrator.u_modified = false
end

function u_modified!(integrator::ODEIntegrator,bool::Bool)
  integrator.u_modified = bool
end

get_proposed_dt(integrator::ODEIntegrator) = integrator.dtpropose
set_proposed_dt!(integrator::ODEIntegrator,dt) = (integrator.dtpropose = dt)

user_cache(integrator::ODEIntegrator) = user_cache(integrator.cache)
u_cache(integrator::ODEIntegrator) = u_cache(integrator.cache)
du_cache(integrator::ODEIntegrator)= du_cache(integrator.cache)
full_cache(integrator::ODEIntegrator) = chain(user_cache(integrator),u_cache(integrator),du_cache(integrator.cache))
default_non_user_cache(integrator::ODEIntegrator) = chain(u_cache(integrator),du_cache(integrator.cache))
function add_tstop!(integrator::ODEIntegrator,t)
  t < integrator.t && error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
  push!(integrator.opts.tstops,t)
end
user_cache(cache::OrdinaryDiffEqCache) = (cache.u,cache.uprev,cache.tmp)

resize!(integrator::ODEIntegrator,i::Int) = resize!(integrator,integrator.cache,i)
function resize!(integrator::ODEIntegrator,cache,i)
  for c in user_cache(integrator)
    resize!(c,i)
  end
  resize_non_user_cache!(integrator,cache,i)
end

resize_non_user_cache!(integrator::ODEIntegrator,i::Int) = resize_non_user_cache!(integrator,integrator.cache,i)
deleteat_non_user_cache!(integrator::ODEIntegrator,i) = deleteat_non_user_cache!(integrator,integrator.cache,i)
addat_non_user_cache!(integrator::ODEIntegrator,i) = addat_non_user_cache!(integrator,integrator.cache,i)

function resize_non_user_cache!(integrator::ODEIntegrator,cache,i)
  for c in default_non_user_cache(integrator)
    resize!(c,i)
  end
end

function resize_non_user_cache!(integrator::ODEIntegrator,cache::Union{Rosenbrock23Cache,Rosenbrock32Cache},i)
  for c in default_non_user_cache(integrator)
    resize!(c,i)
  end
  for c in vecu_cache(integrator.cache)
    resize!(c,i)
  end
  Jvec = vec(cache.J)
  cache.J = reshape(resize!(Jvec,i*i),i,i)
  Wvec = vec(cache.W)
  cache.W = reshape(resize!(Wvec,i*i),i,i)
  resize!(cache.jac_config.duals[1],i)
end
user_cache(cache::Union{Rosenbrock23Cache,Rosenbrock32Cache}) = (cache.u,cache.uprev,cache.jac_config.duals[2])

function resize_non_user_cache!(integrator::ODEIntegrator,cache::Union{GenericImplicitEulerCache,GenericTrapezoidCache},i)
  for c in default_non_user_cache(integrator)
    resize!(c,i)
  end
  for c in vecu_cache(integrator.cache)
    resize!(c,i)
  end
  for c in dual_cache(integrator.cache)
    resize!(c.du,i)
    resize!(c.dual_du,i)
  end
  cache.nl_rhs = integrator.alg.nlsolve(Val{:init},cache.rhs,cache.uhold)
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
  for c in user_cache(integrator)
    deleteat!(c,idxs)
  end
  deleteat_non_user_cache!(integrator,integrator.cache,idxs)
end

function addat!(integrator::ODEIntegrator,idxs)
  for c in user_cache(integrator)
    addat!(c,idxs)
  end
  addat_non_user_cache!(integrator,integrator.cache,idxs)
end

function terminate!(integrator::ODEIntegrator)
  integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
end
