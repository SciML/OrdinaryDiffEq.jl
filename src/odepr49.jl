function solve{uType,tType,isinplace,algType<:ODEIterAlgorithm,F}(prob::AbstractODEProblem{uType,tType,Val{isinplace},F},
    alg::Type{algType},timeseries=[],ts=[],ks=[];dense=true,save_timeseries=true,
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
  if alg <: rk23
    solver = ODE.RKIntegrator{FoA,:rk23}
  elseif alg <: rk45
    solver = ODE.RKIntegrator{FoA,:dopri5}
  elseif alg <: feh78
    solver = ODE.RKIntegrator{FoA,:feh78}
  elseif alg <: ModifiedRosenbrockIntegrator
    solver = ODE.ModifiedRosenbrockIntegrator
  elseif alg <: feuler
    solver = ODE.RKIntegratorFixed{:feuler}
  elseif alg <: midpoint
    solver = ODE.RKIntegratorFixed{:midpoint}
  elseif alg <: heun
    solver = ODE.RKIntegratorFixed{:heun}
  elseif alg <: rk4
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

  saveat_idxs = find((x)->x∈saveat,ts)
  t_nosaveat = view(ts,symdiff(1:length(ts),saveat_idxs))
  u_nosaveat = view(timeseries,symdiff(1:length(ts),saveat_idxs))

  if dense
    interp = (tvals) -> ode_interpolation(tvals,t_nosaveat,u_nosaveat,ks,alg,f!)
  else
    interp = (tvals) -> nothing
  end

  build_ode_solution(prob,alg,ts,timeseries,
                    dense=dense,k=ks,interp=interp,
                    timeseries_errors = timeseries_errors,
                    dense_errors = dense_errors)
end

const ODEJL_OPTION_LIST = Set([:tout,:tstop,:reltol,:abstol,:minstep,:maxstep,:initstep,:norm,:maxiters,:isoutofdomain])
const ODEJL_ALIASES = Dict{Symbol,Symbol}(:minstep=>:dtmin,:maxstep=>:dtmax,:initstep=>:dt,:tstop=>:T,:maxiters=>:maxiters)
const ODEJL_ALIASES_REVERSED = Dict{Symbol,Symbol}([(v,k) for (k,v) in ODEJL_ALIASES])
