"""
`solve(prob::ODEProblem,tspan)`

Solves the ODE defined by prob on the interval tspan. If not given, tspan defaults to [0,1].

Please see the solver documentation.
"""
function solve(prob::AbstractODEProblem,alg=DefaultODEAlgorithm(),timeseries=[],ts=[],ks=[];kwargs...)
  tspan = prob.tspan
  if tspan[end]-tspan[1]<0
    tspan = vec(tspan)
    error("final time must be greater than starting time. Aborting.")
  end
  atomloaded = isdefined(Main,:Atom)
  o = KW(kwargs)
  o[:t] = tspan[1]
  o[:Ts] = tspan[2:end]
  @unpack u0,isinplace = prob
  uType = typeof(u0)
  uEltype = eltype(u0)

  command_opts = copy(DIFFERENTIALEQUATIONSJL_DEFAULT_OPTIONS)
  for (k,v) in o
    command_opts[k]=v
  end
  # Get the control variables
  @unpack save_timeseries, progressbar = command_opts

  if command_opts[:callback] == nothing
    callback = ODE_DEFAULT_CALLBACK
    custom_callback = false
  else
    callback = command_opts[:callback]
    custom_callback = true
  end

  if uEltype<:Number
    u = copy(u0)
  else
    u = deepcopy(u0)
  end

  ks = Vector{uType}(0)

  if typeof(alg) <: OrdinaryDiffEqAlgorithm
    o2 = copy(DIFFERENTIALEQUATIONSJL_DEFAULT_OPTIONS)
    for (k,v) in o
      o2[k]=v
    end
    o = o2

    dt = o[:dt]
    order = alg.order
    adaptiveorder = 0

    if typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm
      adaptiveorder = alg.adaptiveorder
      if o[:adaptive] == true
        dt = 1.0*dt # Convert to float in a way that keeps units
      end
    else
      o[:adaptive] = false
    end

    if typeof(alg) <: ExplicitRK
      @unpack order,adaptiveorder = o[:tableau]
    end

    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du) -> (du[:] = prob.f(t,u))
    else
      f! = prob.f
    end
    if dt==0
      dt = ode_determine_initdt(u0,float(tspan[1]),o[:abstol],o[:reltol],o[:internalnorm],f!,order)
    end

    if o[:tType] == nothing # if tType is not specified, grab it from dt which defaults to 0.0 => Float64
      tType=typeof(dt)
    else
      tType = o[:tType]
    end

    if o[:dtmax] == nothing
      o[:dtmax] = tType((tspan[end]-tspan[1]))
    end
    if o[:dtmin] == nothing
      if tType <: AbstractFloat
        o[:dtmin] = tType(10)*eps(tType)
      else
        o[:dtmin] = tType(1//10^(10))
      end
    end

    if uType <: Number
      uEltypeNoUnits = typeof(u./u)
    else
      uEltypeNoUnits = eltype(u./u)
    end

    Ts = map(tType,o[:Ts])
    t = tType(o[:t])
    rate_prototype = u/zero(t)
    rateType = typeof(rate_prototype) ## Can be different if united

    saveat = tType[convert(tType,x) for x in setdiff(o[:saveat],tspan)]

    if o[:calck]==nothing
      calck = !isempty(saveat) || o[:dense]
    else
      calck = o[:calck]
    end

    ### Algorithm-specific defaults ###

    if o[:qmin] == nothing # Use default qmin
      if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
        qmin = 0.2
      elseif typeof(alg) <: DP8
        qmin = 0.333
      else
        qmin = 0.2
      end
    else
      qmin = o[:qmin]
    end
    if o[:qmax] == nothing # Use default qmax
      if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
        qmax = 10.0
      elseif typeof(alg) <: DP8
        qmax = 6.0
      else
        qmax = 10.0
      end
    else
      qmax = o[:qmax]
    end
    if o[:beta2] == nothing # Use default β₂
      if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
        β₂ = 0.04
      elseif typeof(alg) <: DP8
        β₂ = 0.00
      else
        β₂ = 0.4 / order
      end
    else
      β₂ = o[:beta2]
    end
    if o[:beta1] == nothing # Use default β₁
      if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
        β₁ = 1/order - .75β₂
      elseif typeof(alg) <: DP8
        β₁ = 1/order - .2β₂
      else
        β₁ = .7/order
      end
    else
      β₁ = o[:beta1]
    end
    fsal = false
    if isfsal(alg)
      fsal = true
    elseif typeof(alg) <: ExplicitRK
      @unpack fsal = o[:tableau]
    end

    o[:abstol] = uEltype(1)*o[:abstol]

    if isspecialdense(alg)
      ksEltype = Vector{rateType} # Store more ks for the special algs
    else
      ksEltype = rateType # Makes simple_dense
    end

    timeseries = convert(Vector{uType},timeseries)
    ts = convert(Vector{tType},ts)
    ks = convert(Vector{ksEltype},ks)
    if length(timeseries) == 0
      push!(timeseries,copy(u))
    else
      timeseries[1] = copy(u)
    end

    if length(ts) == 0
      push!(ts,t)
    else
      timeseries[1] = copy(u)
    end

    if ksEltype == rateType
      if uType <: Number
        rate_prototype = f!(t,u)
      else
        f!(t,u,rate_prototype)
      end
      push!(ks,rate_prototype)
    else # Just push a dummy in for special dense since first is not used.
      push!(ks,[rate_prototype])
    end
    γ = o[:gamma]
    dtmax = o[:dtmax]
    dtmin = o[:dtmin]
    @unpack maxiters,timeseries_steps,save_timeseries,adaptive,progress_steps,abstol,reltol,internalnorm,tableau,autodiff,qoldinit,dense = o
    # @code_warntype ode_solve(ODEIntegrator{alg,uType,uEltype,ndims(u)+1,tType,uEltypeNoUnits,rateType,ksEltype}(timeseries,ts,ks,f!,u,t,k,dt,Ts,maxiters,timeseries_steps,save_timeseries,adaptive,abstol,reltol,γ,qmax,qmin,dtmax,dtmin,internalnorm,progressbar,tableau,autodiff,adaptiveorder,order,atomloaded,progress_steps,β₁,β₂,qoldinit,fsal,dense,saveat,alg,callback,custom_callback,calck))
    u,t = ode_solve(ODEIntegrator{typeof(alg),uType,uEltype,ndims(u)+1,tType,uEltypeNoUnits,rateType,ksEltype}(timeseries,ts,ks,f!,u,t,dt,Ts,maxiters,timeseries_steps,save_timeseries,adaptive,abstol,reltol,γ,qmax,qmin,dtmax,dtmin,internalnorm,progressbar,tableau,autodiff,adaptiveorder,order,atomloaded,progress_steps,β₁,β₂,qoldinit,fsal,dense,saveat,alg,callback,custom_callback,calck))
    if ts[end] != t
      push!(ts,t)
      push!(timeseries,u)
    end
  elseif typeof(alg) <: ODEInterfaceAlgorithm
    sizeu = size(u)
    if typeof(u) <: Number
      u = [u]
    end
    o[:Ts] = float(o[:Ts])
    o[:t] = float(o[:t])
    t = o[:t]; Ts = o[:Ts]
    saveat = [float(x) for x in command_opts[:saveat]]
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du) -> (du[:] = vec(prob.f(t,reshape(u,sizeu))); nothing)
    else
      f! = (t,u,du) -> (prob.f(t,reshape(u,sizeu),reshape(du,sizeu)); u = vec(u); du=vec(du); nothing)
    end
    initialize_backend(:ODEInterface)
    o[:RHS_CALLMODE] = ODEInterface.RHS_CALL_INSITU
    dict = buildOptions(o,ODEINTERFACE_OPTION_LIST,ODEINTERFACE_ALIASES,ODEINTERFACE_ALIASES_REVERSED)
    opts = ODEInterface.OptionsODE([Pair(ODEINTERFACE_STRINGS[k],v) for (k,v) in dict]...) #Convert to the strings
    du = similar(u)
    if typeof(alg) <: dopri5
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.dopri5,f!,[t;Ts],vec(u),opts)
    elseif typeof(alg) <: dop853
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.dop853,f!,[t;Ts],vec(u),opts)
    elseif typeof(alg) <: odex
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.odex,f!,[t;Ts],vec(u),opts)
    elseif typeof(alg) <: seulex
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.seulex,f!,[t;Ts],vec(u),opts)
    elseif typeof(alg) <: radau
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.radau,f!,[t;Ts],vec(u),opts)
    elseif typeof(alg) <: radau5
      ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.radau5,f!,[t;Ts],vec(u),opts)
    end
    if retcode < 0
      if retcode == -1
        warn("Input is not consistent.")
      elseif retcode == -2
        warn("Interrupted. Larger maxiters is needed.")
      elseif retcode == -3
        warn("Step size went too small.")
      elseif retcode == -4
        warn("Interrupted. Problem is probably stiff.")
      end
    end
    t = ts[end]
    if typeof(u0)<:AbstractArray
      timeseries = Vector{uType}(0)
      for i=1:size(vectimeseries,1)
        push!(timeseries,reshape(view(vectimeseries,i,:,)',sizeu))
      end
    else
      timeseries = vec(vectimeseries)
    end
    u = timeseries[end]
  elseif typeof(alg) <: ODEIterAlgorithm
    if typeof(u) <: Number
      u = [u]
    end
    # Needs robustness
    o[:Ts] = float(o[:Ts])
    o[:t] = float(o[:t])
    t = o[:t]; Ts = o[:Ts]
    o[:T] = Ts[end]
    saveat = [float(x) for x in command_opts[:saveat]]
    initialize_backend(:ODEJL)
    opts = buildOptions(o,ODEJL_OPTION_LIST,ODEJL_ALIASES,ODEJL_ALIASES_REVERSED)
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du) -> (du[:] = prob.f(t,u))
    else
      f! = prob.f
    end
    ode  = ODE.ExplicitODE(t,u,f!)
    # adaptive==true ? FoA=:adaptive : FoA=:fixed #Currently limied to only adaptive
    FoA = :adaptive
    if typeof(alg) <: ode23
      solver = ODE.RKIntegrator{FoA,:rk23}
    elseif typeof(alg) <: ode45
      solver = ODE.RKIntegrator{FoA,:dopri5}
    elseif typeof(alg) <: ode78
      solver = ODE.RKIntegrator{FoA,:feh78}
    elseif typeof(alg) <: ode23s
      solver = ODE.ModifiedRosenbrockIntegrator
    elseif typeof(alg) <: ode1
      solver = ODE.RKIntegratorFixed{:feuler}
    elseif typeof(alg) <: ode2_midpoint
      solver = ODE.RKIntegratorFixed{:midpoint}
    elseif typeof(alg) <: ode2_heun
      solver = ODE.RKIntegratorFixed{:heun}
    elseif typeof(alg) <: ode4
      solver = ODE.RKIntegratorFixed{:rk4}
    elseif typeof(alg) <: ode45_fe
      solver = ODE.RKIntegrator{FoA,:rk45}
    end
    out = ODE.solve(ode;solver=solver,opts...)
    timeseries = out.y
    ts = out.t
    ks = out.dy
    if length(out.y[1])==1
      tmp = Vector{eltype(out.y[1])}(length(out.y))
      tmp_dy = Vector{eltype(out.dy[1])}(length(out.dy))
      for i in 1:length(out.y)
        tmp[i] = out.y[i][1]
        tmp_dy[i] = out.dy[i][1]
      end
      timeseries = tmp
      ks = tmp_dy
    end
    t = ts[end]
    u = timeseries[end]
  elseif typeof(alg) <: SundialsAlgorithm
    if typeof(alg) == cvode_BDF
      integrator = :BDF
    elseif typeof(alg) ==  cvode_Adams
      integrator = :Adams
    end

    sizeu = size(u)
    if typeof(u) <: Number
      u = [u]
    end
    u = map(Float64,u) # Needs Float64
    # Needs robustness
    o[:Ts] = map(Float64,o[:Ts])
    o[:t] = map(Float64,o[:t])
    t = o[:t]; Ts = o[:Ts];
    saveat = [float(x) for x in command_opts[:saveat]]
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du) -> (du[:] = vec(prob.f(t,reshape(u,sizeu))); 0)
    else
      f! = (t,u,du) -> (prob.f(t,reshape(u,sizeu),reshape(du,sizeu)); u = vec(u); du=vec(du); 0)
    end
    ts = [t;Ts]
    @unpack abstol, reltol = command_opts
    if command_opts[:adaptive]
      ts, vectimeseries = Sundials.cvode_fulloutput(f!,vec(u),ts;integrator=integrator,abstol=float(abstol),reltol=float(reltol))
      timeseries = Vector{uType}(0)
      if typeof(u0)<:AbstractArray
        for i=1:size(vectimeseries,1)
          push!(timeseries,reshape(vectimeseries[i],sizeu))
        end
      else
        for i=1:size(vectimeseries,1)
          push!(timeseries,vectimeseries[i][1])
        end
      end
    else
      dt = command_opts[:dt]
      ts = float(collect(t:dt:Ts[end]))
      if length(Ts)>1
        ts = float([ts;Ts[1:end-1]])
        sort(ts)
      end
      vectimeseries = Sundials.cvode(f!,vec(u),ts,integrator=integrator,abstol=float(abstol),reltol=float(reltol))
      timeseries = Vector{uType}(0)
      if typeof(u0)<:AbstractArray
        for i=1:size(vectimeseries,1)
          push!(timeseries,reshape(view(vectimeseries,i,:),sizeu))
        end
      else
        for i=1:size(vectimeseries,1)
          push!(timeseries,vectimeseries[i])
        end
      end
    end
    t = ts[end]
    u = timeseries[end]
  end

  (atomloaded && progressbar) ? Main.Atom.progress(1) : nothing #Use Atom's progressbar if loaded

  if typeof(prob) <: ODETestProblem
    u_analytic = prob.analytic(t,u0)
    timeseries_analytic = Vector{uType}(0)
    for i in 1:size(timeseries,1)
      push!(timeseries_analytic,prob.analytic(ts[i],u0))
    end
    return(ODESolution(u,u_analytic,prob,alg,timeseries=timeseries,t=ts,timeseries_analytic=timeseries_analytic,k=ks,saveat=saveat,
    timeseries_errors = command_opts[:timeseries_errors],
    dense_errors = command_opts[:dense_errors]))
  else
    return(ODESolution(u,prob,alg,timeseries=timeseries,t=ts,k=ks,saveat=saveat))
  end
end

function buildOptions(o,optionlist,aliases,aliases_reversed)
  dict1 = Dict{Symbol,Any}([Pair(k,o[k]) for k in (keys(o) ∩ optionlist)])
  dict2 = Dict([Pair(aliases_reversed[k],o[k]) for k in (keys(o) ∩ values(aliases))])
  merge(dict1,dict2)
end

function ode_determine_initdt(u0,t,abstol,reltol,internalnorm,f,order)
  f₀ = similar(u0); f₁ = similar(u0); u₁ = similar(u0)
  d₀ = internalnorm(u0./(abstol+u0*reltol))
  f(t,u0,f₀)
  d₁ = internalnorm(f₀./(abstol+u0*reltol))
  if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    dt₀ = 1//10^(6)
  else
    dt₀ = (d₀/d₁)/100
  end
  @inbounds for i in eachindex(u0)
     u₁[i] = u0[i] + dt₀*f₀[i]
  end
  f(t+dt₀,u₁,f₁)
  d₂ = internalnorm((f₁.-f₀)./(abstol+u0*reltol))/dt₀
  if max(d₁,d₂)<=1//10^(15)
    dt₁ = max(1//10^(6),dt₀*1//10^(3))
  else
    dt₁ = 10.0^(-(2+log10(max(d₁,d₂)))/(order))
  end
  dt = min(100dt₀,dt₁)
end

function ode_determine_initdt(u0::Number,t,abstol,reltol,internalnorm,f,order)
  d₀ = abs(u0./(abstol+u0*reltol))
  f₀ =f(t,u0)
  d₁ = abs(f₀./(abstol+u0*reltol))
  if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    dt₀ = 1//10^(6)
  else
    dt₀ = (d₀/d₁)/100
  end
  u₁ = u0 + dt₀*f₀
  f₁ = f(t+dt₀,u₁)
  d₂ = abs((f₁-f₀)./(abstol+u0*reltol))/dt₀
  if max(d₁,d₂)<=1//10^(15)
    dt₁ = max(1//10^(6),dt₀*1//10^(3))
  else
    dt₁ = 10.0^(-(2+log10(max(d₁,d₂)))/(order))
  end
  dt = min(100dt₀,dt₁)
end

function plan_ode(alg_hint,abstol,reltol)
  :DP5
end
