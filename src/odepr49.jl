function solve{uType,tType,isinplace,algType<:ODEIterAlgorithm,F}(prob::AbstractODEProblem{uType,tType,Val{isinplace},F},
    alg::Type{algType},timeseries=[],ts=[],ks=[];
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

const ODEJL_OPTION_LIST = Set([:tout,:tstop,:reltol,:abstol,:minstep,:maxstep,:initstep,:norm,:maxiters,:isoutofdomain])
const ODEJL_ALIASES = Dict{Symbol,Symbol}(:minstep=>:dtmin,:maxstep=>:dtmax,:initstep=>:dt,:tstop=>:T,:maxiters=>:maxiters)
const ODEJL_ALIASES_REVERSED = Dict{Symbol,Symbol}([(v,k) for (k,v) in ODEJL_ALIASES])
