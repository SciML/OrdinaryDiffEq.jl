function solve{uType,tType,isinplace,T<:OrdinaryDiffEqAlgorithm}(prob::AbstractODEProblem{uType,tType,
  Val{isinplace}},algType::Type{T}=DefaultODEAlgorithm(),
  timeseries=[],ts=[],ks=[];dt = 0.0,save_timeseries = true,
  timeseries_steps = 1,tableau = ODE_DEFAULT_TABLEAU,
  dense = true,calck = nothing,alg_hint = :nonstiff,
  timeseries_errors = true,dense_errors = false,saveat = Float64[],
  adaptive = true,gamma=.9,abstol=1//10^6,reltol=1//10^3,
  qmax=nothing,qmin=nothing,qoldinit=1//10^4, fullnormalize=true,
  beta2=nothing,beta1=nothing,maxiters = 10000,
  dtmax=tType((prob.tspan[end]-prob.tspan[1])),
  dtmin=tType <: AbstractFloat ? tType(10)*eps(tType) : tType(1//10^(10)),
  autodiff=true,internalnorm = ODE_DEFAULT_NORM,
  progressbar=false,progress_steps=1000,callback=nothing,kwargs...)

  alg = algType()

  tspan = prob.tspan

  if tspan[end]-tspan[1]<tType(0)
    error("final time must be greater than starting time. Aborting.")
  end
  atomloaded = isdefined(Main,:Atom)

  t = tspan[1]
  Ts = tspan[2:end]
  u0 = prob.u0
  uEltype = eltype(u0)

  # Get the control variables

  if callback == nothing
    callback = ODE_DEFAULT_CALLBACK
    custom_callback = false
  else
    custom_callback = true
  end

  if uEltype<:Number
    u = copy(u0)
  else
    u = deepcopy(u0)
  end

  ks = Vector{uType}(0)

  order = alg.order
  adaptiveorder = 0

  if typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm
    adaptiveorder = alg.adaptiveorder
    if adaptive == true
      dt = 1.0*dt # Convert to float in a way that keeps units
    end
  else
    adaptive = false
  end

  if typeof(alg) <: ExplicitRK
    @unpack order,adaptiveorder = tableau
  end

  if !isinplace && typeof(u)<:AbstractArray
    f! = (t,u,du) -> (du[:] = prob.f(t,u))
  else
    f! = prob.f
  end

  if uType <: Number
    uEltypeNoUnits = typeof(u./u)
  else
    uEltypeNoUnits = eltype(u./u)
  end

  if dt==0
    dt = ode_determine_initdt(u0,t,uEltype(abstol),uEltypeNoUnits(reltol),internalnorm,f!,order)
  end

  rate_prototype = u/zero(t)
  rateType = typeof(rate_prototype) ## Can be different if united

  saveat = tType[convert(tType,x) for x in setdiff(saveat,tspan)]

  if calck==nothing
    calck = !isempty(saveat) || dense
  end

  ### Algorithm-specific defaults ###

  if qmin == nothing # Use default qmin
    if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
      qmin = 0.2
    elseif typeof(alg) <: DP8
      qmin = 0.333
    else
      qmin = 0.2
    end
  end
  if qmax == nothing # Use default qmax
    if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
      qmax = 10.0
    elseif typeof(alg) <: DP8
      qmax = 6.0
    else
      qmax = 10.0
    end
  end
  if beta2 == nothing # Use default β₂
    if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
      β₂ = 0.04
    elseif typeof(alg) <: DP8
      β₂ = 0.00
    else
      β₂ = 0.4 / order
    end
  else
    β₂ = beta2
  end
  if beta1 == nothing # Use default β₁
    if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded
      β₁ = 1/order - .75β₂
    elseif typeof(alg) <: DP8
      β₁ = 1/order - .2β₂
    else
      β₁ = .7/order
    end
  else
    β₁ = beta1
  end

  fsal = false
  if isfsal(alg)
    fsal = true
  elseif typeof(alg) <: ExplicitRK
    @unpack fsal = tableau
  end

  abstol = uEltype(1)*abstol

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
  γ = gamma
  # @code_warntype ode_solve(ODEIntegrator{alg,uType,uEltype,ndims(u)+1,tType,uEltypeNoUnits,rateType,ksEltype}(timeseries,ts,ks,f!,u,t,k,dt,Ts,maxiters,timeseries_steps,save_timeseries,adaptive,abstol,reltol,γ,qmax,qmin,dtmax,dtmin,internalnorm,progressbar,tableau,autodiff,adaptiveorder,order,atomloaded,progress_steps,β₁,β₂,qoldinit,fsal,dense,saveat,alg,callback,custom_callback,calck))
  u,t = ode_solve(ODEIntegrator{typeof(alg),uType,uEltype,ndims(u)+1,tType,uEltypeNoUnits,rateType,ksEltype}(timeseries,ts,ks,f!,u,t,dt,Ts,maxiters,timeseries_steps,save_timeseries,adaptive,abstol,reltol,γ,qmax,qmin,dtmax,dtmin,internalnorm,progressbar,tableau,autodiff,adaptiveorder,order,atomloaded,progress_steps,β₁,β₂,qoldinit,fsal,dense,saveat,alg,callback,custom_callback,calck))

  if typeof(prob) <: ODETestProblem
    timeseries_analytic = Vector{uType}(0)
    for i in 1:size(timeseries,1)
      push!(timeseries_analytic,prob.analytic(ts[i],u0))
    end
    return(ODESolution(ts,timeseries,prob,alg,
    u_analytic=timeseries_analytic,
    k=ks,saveat=saveat,
    timeseries_errors = timeseries_errors,
    dense_errors = dense_errors))
  else
    return(ODESolution(ts,timeseries,prob,alg,k=ks,saveat=saveat))
  end
end



function solve{uType,tType,isinplace,T<:ODEInterfaceAlgorithm}(
    prob::AbstractODEProblem{uType,tType,Val{isinplace}},
    alg::Type{T},timeseries=[],ts=[],ks=[];kwargs...)

  tspan = prob.tspan

  if tspan[end]-tspan[1]<tType(0)
    error("final time must be greater than starting time. Aborting.")
  end

  o = KW(kwargs)

  u0 = prob.u0

  if typeof(u0) <: Number
    u = [u0]
  else
    u = deepcopy(u0)
  end

  sizeu = size(u)

  if !isinplace && typeof(u)<:AbstractArray
    f! = (t,u,du) -> (du[:] = vec(prob.f(t,reshape(u,sizeu))); nothing)
  else
    f! = (t,u,du) -> (prob.f(t,reshape(u,sizeu),reshape(du,sizeu)); u = vec(u); du=vec(du); nothing)
  end

  initialize_backend(:ODEInterface)
  o[:RHS_CALLMODE] = ODEInterface.RHS_CALL_INSITU
  dict = buildOptions(o,ODEINTERFACE_OPTION_LIST,ODEINTERFACE_ALIASES,ODEINTERFACE_ALIASES_REVERSED)
  opts = ODEInterface.OptionsODE([Pair(ODEINTERFACE_STRINGS[k],v) for (k,v) in dict]...) #Convert to the strings
  if alg <: dopri5
    ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.dopri5,f!,tspan,vec(u),opts)
  elseif alg <: dop853
    ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.dop853,f!,tspan,vec(u),opts)
  elseif alg <: odex
    ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.odex,f!,tspan,vec(u),opts)
  elseif alg <: seulex
    ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.seulex,f!,tspan,vec(u),opts)
  elseif alg <: radau
    ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.radau,f!,tspan,vec(u),opts)
  elseif alg <: radau5
    ts,vectimeseries,retcode,stats = ODEInterface.odecall(ODEInterface.radau5,f!,tspan,vec(u),opts)
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

  if typeof(u0)<:AbstractArray
    timeseries = Vector{uType}(0)
    for i=1:size(vectimeseries,1)
      push!(timeseries,reshape(view(vectimeseries,i,:,)',sizeu))
    end
  else
    timeseries = vec(vectimeseries)
  end

  if typeof(prob) <: ODETestProblem
    timeseries_analytic = Vector{uType}(0)
    for i in 1:size(timeseries,1)
      push!(timeseries_analytic,prob.analytic(ts[i],u0))
    end
    return(ODESolution(ts,timeseries,prob,alg,u_analytic=timeseries_analytic))
  else
    return(ODESolution(ts,timeseries,prob,alg))
  end
end

function solve{uType,tType,isinplace,algType<:SundialsAlgorithm}(prob::AbstractODEProblem{uType,tType,Val{isinplace}},
    alg::Type{algType}=DefaultODEAlgorithm(),timeseries=[],ts=[],ks=[];
    callback=()->nothing,abstol=1/10^6,reltol=1/10^3,saveat=Float64[],adaptive=true,
    timeseries_errors=true,dense_errors=false,save_timeseries=true,
    kwargs...)
  tspan = prob.tspan

  atomloaded = isdefined(Main,:Atom)
  u0 = prob.u0

  if typeof(u0) <: Number
    u = [u0]
  else
    u = deepcopy(u0)
  end
  if alg == cvode_BDF
    integrator = :BDF
  elseif alg ==  cvode_Adams
    integrator = :Adams
  end

  sizeu = size(u)
  if !isinplace && typeof(u)<:AbstractArray
    f! = (t,u,du) -> (du[:] = vec(prob.f(t,reshape(u,sizeu))); 0)
  else
    f! = (t,u,du) -> (prob.f(t,reshape(u,sizeu),reshape(du,sizeu)); u = vec(u); du=vec(du); 0)
  end

  ts = sort([tspan;saveat])
  sort(ts)

  if save_timeseries
    ts, vectimeseries = Sundials.cvode_fulloutput(f!,vec(u),ts;integrator=integrator,abstol=abstol,reltol=reltol)
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
    vectimeseries = Sundials.cvode(f!,vec(u),ts,integrator=integrator,abstol=abstol,reltol=reltol)
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

  if typeof(prob) <: ODETestProblem
    timeseries_analytic = Vector{uType}(0)
    for i in 1:size(timeseries,1)
      push!(timeseries_analytic,prob.analytic(ts[i],u0))
    end
    return(ODESolution(ts,timeseries,prob,alg,
    u_analytic=timeseries_analytic,
    timeseries_errors = timeseries_errors,
    dense_errors = dense_errors))
  else
    return(ODESolution(ts,timeseries,prob,alg))
  end
end

function buildOptions(o,optionlist,aliases,aliases_reversed)
  dict1 = Dict{Symbol,Any}([Pair(k,o[k]) for k in (keys(o) ∩ optionlist)])
  dict2 = Dict([Pair(aliases_reversed[k],o[k]) for k in (keys(o) ∩ values(aliases))])
  merge(dict1,dict2)
end

function ode_determine_initdt{uType,tType,uEltypeNoUnits}(u0::uType,t::tType,abstol,reltol::uEltypeNoUnits,internalnorm,f,order)
  f₀ = similar(u0./t); f₁ = similar(u0./t); u₁ = similar(u0)
  d₀ = internalnorm(u0./(abstol+u0*reltol))
  f(t,u0,f₀)
  d₁ = internalnorm(f₀./(abstol+u0*reltol)*tType(1))/tType(1)
  T0 = typeof(d₀)
  T1 = typeof(d₁)
  if d₀ < T0(1//10^(5)) || d₁ < T1(1//10^(5))
    dt₀ = tType(1//10^(6))
  else
    dt₀ = tType((d₀/d₁)/100)
  end
  @inbounds for i in eachindex(u0)
     u₁[i] = u0[i] + dt₀*f₀[i]
  end
  f(t+dt₀,u₁,f₁)
  d₂ = internalnorm((f₁.-f₀)./(abstol+u0*reltol)*tType(1))/dt₀
  if max(d₁,d₂)<=T1(1//10^(15))
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(max(d₁,d₂)/T1(1)))/(order)))
  end
  dt = min(100dt₀,dt₁)
end

function ode_determine_initdt{uType<:Number,tType,uEltypeNoUnits}(u0::uType,t::tType,abstol,reltol::uEltypeNoUnits,internalnorm,f,order)
  d₀ = abs(u0./(abstol+u0*reltol))
  f₀ = f(t,u0)
  d₁ = abs(f₀./(abstol+u0*reltol))
  T0 = typeof(d₀)
  T1 = typeof(d₁)
  if d₀ < T0(1//10^(5)) || d₁ < T1(1//10^(5))
    dt₀ = tType(1//10^(6))
  else
    dt₀ = tType((d₀/d₁)/100)
  end
  u₁ = u0 + dt₀*f₀
  f₁ = f(t+dt₀,u₁)
  d₂ = abs((f₁-f₀)./(abstol+u0*reltol))/dt₀*tType(1)
  if max(d₁,d₂) <= T1(1//10^(15))
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(max(d₁,d₂)/T1(1)))/(order)))
  end
  dt = min(100dt₀,dt₁)
end




function solve{uType,tType,isinplace,algType<:ODEIterAlgorithm}(prob::AbstractODEProblem{uType,tType,Val{isinplace}},
    alg::Type{algType}=DefaultODEAlgorithm(),timeseries=[],ts=[],ks=[];
    saveat=[],callback=()->nothing,timeseries_errors=true,dense_errors=false,
    kwargs...)
  tspan = prob.tspan

  if tspan[end]-tspan[1]<tType(0)
    error("final time must be greater than starting time. Aborting.")
  end
  atomloaded = isdefined(Main,:Atom)
  o = KW(kwargs)
  t = tspan[1]
  u0 = prob.u0
  o[:T] = tspan[end]

  if typeof(u0) <: Number
    u = [u0]
  else
    u = deepcopy(u0)
  end

  sizeu = size(u)

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
  if alg <: ode23
    solver = ODE.RKIntegrator{FoA,:rk23}
  elseif alg <: rk45
    solver = ODE.RKIntegrator{FoA,:dopri5}
  elseif alg <: ode78
    solver = ODE.RKIntegrator{FoA,:feh78}
  elseif alg <: ode23s
    solver = ODE.ModifiedRosenbrockIntegrator
  elseif alg <: ode1
    solver = ODE.RKIntegratorFixed{:feuler}
  elseif alg <: ode2_midpoint
    solver = ODE.RKIntegratorFixed{:midpoint}
  elseif alg <: ode2_heun
    solver = ODE.RKIntegratorFixed{:heun}
  elseif alg <: ode4
    solver = ODE.RKIntegratorFixed{:rk4}
  elseif alg <: feh45
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

  if typeof(prob) <: ODETestProblem
    timeseries_analytic = Vector{uType}(0)
    for i in 1:size(timeseries,1)
      push!(timeseries_analytic,prob.analytic(ts[i],u0))
    end
    return(ODESolution(ts,timeseries,prob,alg,
    u_analytic=timeseries_analytic,
    saveat=saveat,
    timeseries_errors = timeseries_errors,
    dense_errors = dense_errors))
  else
    return(ODESolution(ts,timeseries,prob,alg,saveat=saveat))
  end
end
