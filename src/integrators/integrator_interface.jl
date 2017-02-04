@inline function change_t_via_interpolation!{T}(integrator,t,modify_save_endpoint::Type{Val{T}}=Val{false})
  # Can get rid of an allocation here with a function
  # get_tmp_arr(integrator.cache) which gives a pointer to some
  # cache array which can be modified.
  if t < integrator.tprev
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

@inline function reeval_internals_due_to_modification!(integrator)
  if integrator.opts.calck
    resize!(integrator.k,integrator.kshortsize) # Reset k for next step!
    ode_addsteps!(integrator,integrator.f,Val{true},Val{false})
  end
  integrator.u_modified = false
end

@inline function u_modified!(integrator::ODEIntegrator,bool::Bool)
  integrator.u_modified = bool
end

u_cache(integrator::ODEIntegrator) = u_cache(integrator.cache)
du_cache(integrator::ODEIntegrator)= du_cache(integrator.cache)
full_cache(integrator::ODEIntegrator) = chain(u_cache(integrator),du_cache(integrator.cache))
add_tstop!(integrator::ODEIntegrator,t) = push!(integrator.opts.tstops,t)
