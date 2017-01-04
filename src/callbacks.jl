@inline function determine_event_occurance(integrator,callback)
  event_occurred = false
  if callback.interp_points!=0
    ode_addsteps!(integrator)
  end
  Θs = linspace(typeof(integrator.t)(0),typeof(integrator.t)(1),callback.interp_points)
  interp_index = 0
  # Check if the event occured
  if callback.condition(integrator.t,integrator.u,integrator)<=0
    event_occurred = true
    interp_index = callback.interp_points
  elseif callback.interp_points!=0 # Use the interpolants for safety checking
    for i in 2:length(Θs)-1
      if callback.condition(integrator.t+integrator.dt*Θs[i],ode_interpolant(Θs[i],integrator),integrator)<0
        event_occurred = true
        interp_index = i
        break
      end
    end
  end
  event_occurred,interp_index,Θs
end

function apply_callback!(integrator,callback)
  event_occurred,interp_index,Θs = determine_event_occurance(integrator,callback)
  if event_occurred
    if callback.interp_points!=0
      top_Θ = Θs[interp_index] # Top at the smallest
    else
      top_Θ = typeof(integrator.t)(1)
    end
    if callback.rootfind
      find_zero = (Θ) -> begin
        callback.condition(integrator.tprev+Θ*integrator.dt,ode_interpolant(Θ,integrator),integrator)
      end
      res = prevfloat(fzero(find_zero,typeof(integrator.t)(0),top_Θ))
      val = ode_interpolant(res,integrator)
      recursivecopy!(integrator.u,val)
      integrator.dt *= res
    elseif interp_index != callback.interp_points
        integrator.dt *= Θs[interp_index]
        recursivecopy!(integrator.u,ode_interpolant(Θs[interp_index],integrator))
    end
    # If no solve and no interpolants, just use endpoint

    integrator.t = integrator.tprev + integrator.dt
    if integrator.opts.calck
      if isspecialdense(integrator.alg)
        resize!(integrator.k,integrator.kshortsize) # Reset k for next step!
        ode_addsteps!(integrator,Val{true},Val{false})
      elseif typeof(integrator.u) <: Number
        integrator.k = integrator.f(integrator.t,integrator.u)
      else
        integrator.f(integrator.t,integrator.u,integrator.k)
      end
    end
  end

  if callback.save_positions[1]
    savevalues!(integrator)
  end

  if event_occurred
    callback.affect!(integrator)
    if integrator.opts.calck
      if isspecialdense(integrator.alg)
        resize!(integrator.k,integrator.kshortsize) # Reset k for next step!
        ode_addsteps!(integrator,Val{true},Val{false})
      else
        if typeof(integrator.u) <: Number
          integrator.k = integrator.f(integrator.t,integrator.u)
        else
          integrator.f(integrator.t,integrator.u,integrator.k)
        end
      end
    end

    if callback.save_positions[2]
      savevalues!(integrator)
    end

    if (isfsal(integrator.alg) && !isspecialdense(integrator.alg)) || (isspecialdense(integrator.alg) && isfsal(integrator.alg) && !(integrator.fsalfirst===integrator.k)) ## This will stop double compute for simpledense FSAL
      integrator.reeval_fsal = true
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
