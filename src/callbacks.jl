macro ode_callback(ex)
  esc(quote
    function (t,cache,T,Ts,integrator)
      event_occurred = false
      $(ex)
      t,T
    end
  end)
end

macro ode_event(event_f,apply_event!,rootfind_event_loc=true,interp_points=5,terminate_on_event=false,dt_safety=1)
  esc(quote
    # Event Handling
    if $interp_points!=0
      ode_addsteps!(integrator.k,integrator.tprev,integrator.uprev,integrator.dt,integrator.alg,integrator.f)
      Θs = linspace(typeof(t)(0),typeof(t)(1),$(interp_points))
    end
    interp_index = 0
    # Check if the event occured
    if $event_f(integrator.t,integrator.u)<=0
      event_occurred = true
      interp_index = $interp_points
    elseif $interp_points!=0 # Use the interpolants for safety checking
      for i in 2:length(Θs)-1
        if $event_f(integrator.t+integrator.dt*Θs[i],ode_interpolant(Θs[i],integrator.dt,integrator.uprev,integrator.u,integrator.kprev,integrator.k,integrator.alg))<0
          event_occurred = true
          interp_index = i
          break
        end
      end
    end

    if event_occurred
      top_Θ = Θs[interp_index] # Top at the smallest
      if $rootfind_event_loc
        find_zero = (Θ) -> begin
          $event_f(integrator.tprev+Θ*integrator.dt,ode_interpolant(Θ,integrator.dt,integrator.uprev,integrator.u,integrator.kprev,integrator.k,integrator.alg))
        end
        res = prevfloat(fzero(find_zero,typeof(t)(0),top_Θ))
        val = ode_interpolant(res,integrator.dt,integrator.uprev,integrator.u,integrator.kprev,integrator.k,integrator.alg)
        copy!(integrator.u,val)
        integrator.dt *= res
      elseif interp_index != $interp_points
          integrator.dt *= Θs[interp_index]
          copy!(integrator.u,ode_interpolant(Θs[interp_index],integrator.dt,integrator.uprev,integrator.u,integrator.kprev,integrator.k,integrator.alg))
      end
      # If no solve and no interpolants, just use endpoint

      integrator.t = integrator.tprev + integrator.dt

      if integrator.opts.calck
        if isspecialdense(integrator.alg)
          resize!(integrator.k,integrator.kshortsize) # Reset k for next step
          integrator.k = typeof(integrator.k)() # Make a local blank k for saving
          ode_addsteps!(integrator.k,integrator.tprev,integrator.uprev,integrator.dt,integrator.alg,integrator.f)
        elseif typeof(integrator.u) <: Number
          integrator.k = integrator.f(integrator.t,integrator.u)
        else
          integrator.f(integrator.t,integrator.u,integrator.k)
        end
      end
    end

    ode_savevalues!(integrator)
    if event_occurred
      if $terminate_on_event
        @ode_terminate
      else
        $apply_event!(integrator.u,cache)
        if integrator.opts.calck
          if !isspecialdense(integrator.alg)
            if typeof(integrator.u) <: Number
              integrator.k = integrator.f(integrator.t,integrator.u)
            else
              integrator.f(integrator.t,integrator.u,integrator.k)
            end
          end
        end
        ode_savevalues!(integrator)
        if integrator.fsal
          integrator.reeval_fsal = true
        end
        integrator.dt_mod *= $dt_safety # Safety dt change
      end
    end
  end)
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

@def ode_terminate begin
  T = integrator.t
  while length(Ts)>1
    pop!(Ts)
  end
end
