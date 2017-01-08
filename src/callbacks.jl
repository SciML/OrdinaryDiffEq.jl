@inline function determine_event_occurance(integrator,callback)
  event_occurred = false
  if callback.interp_points!=0
    ode_addsteps!(integrator)
  end
  Θs = linspace(typeof(integrator.t)(0),typeof(integrator.t)(1),callback.interp_points)
  interp_index = 0
  # Check if the event occured
  if typeof(callback.condition) <: Void
    event_occurred = true
  else
    previous_condition = callback.condition(integrator.tprev,integrator.uprev,integrator)
    if isapprox(previous_condition,0,rtol=callback.reltol,atol=callback.abstol)
      prev_sign = 0
    else
      prev_sign = sign(previous_condition)
    end
    if prev_sign*sign(callback.condition(integrator.t,integrator.u,integrator))<0
      event_occurred = true
      interp_index = callback.interp_points
    elseif callback.interp_points!=0 # Use the interpolants for safety checking
      for i in 2:length(Θs)-1
        if prev_sign*callback.condition(integrator.t+integrator.dt*Θs[i],ode_interpolant(Θs[i],integrator),integrator)<0
          event_occurred = true
          interp_index = i
          break
        end
      end
    end
  end
  event_occurred,interp_index,Θs
end

function apply_callback!(integrator,callback)
  event_occurred,interp_index,Θs = determine_event_occurance(integrator,callback)
  if event_occurred
    if typeof(callback.condition) <: Void
      new_t = integrator.t
    else
      if callback.interp_points!=0
        top_Θ = Θs[interp_index] # Top at the smallest
      else
        top_Θ = typeof(integrator.t)(1)
      end
      if callback.rootfind
        find_zero = (Θ) -> begin
          callback.condition(integrator.tprev+Θ*integrator.dt,ode_interpolant(Θ,integrator),integrator)
        end
        Θ = prevfloat(prevfloat(fzero(find_zero,typeof(integrator.t)(0),top_Θ)))
        new_t = integrator.tprev + integrator.dt*Θ
      elseif interp_index != callback.interp_points
        new_t = integrator.tprev + integrator.dt*Θs[interp_index]
      else
        # If no solve and no interpolants, just use endpoint
        new_t = integrator.t
      end
    end
    change_t_via_interpolation!(integrator,new_t)
  end

  if callback.save_positions[1]
    savevalues!(integrator)
  end

  if event_occurred
    integrator.u_modified = true
    callback.affect!(integrator)
    if integrator.u_modified
      reeval_internals_due_to_modification!(integrator)
    end
    if callback.save_positions[2]
      savevalues!(integrator)
    end
  end
end

macro ode_change_cachesize(cache,resize_ex)
  resize_ex = cache_replace_length(resize_ex)
  esc(quote
    for i in 1:length($cache)
      resize!($cache[i],$resize_ex)
    end
  end)
end

macro ode_change_deleteat(cache,deleteat_ex)
  deleteat_ex = cache_replace_length(deleteat_ex)
  esc(quote
    for i in 1:length($cache)
      deleteat!($cache[i],$deleteat_ex)
    end
  end)
end

function cache_replace_length(ex::Expr)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      cache_replace_length(ex)
    elseif isa(arg,Symbol)
      if arg == :length
        ex.args[i] = :(length(cache[i]))
      end
    end
  end
  ex
end

function cache_replace_length(ex::Symbol)
  if ex == :length
    ex = :(length(cache[i]))
  end
  ex
end

function cache_replace_length(ex::Any)
  ex
end

function terminate!(integrator::ODEIntegrator)
  integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
end
