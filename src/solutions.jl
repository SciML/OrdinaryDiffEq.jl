"""
`ODESolution`

Holds the data for the solution to an ODE problem.

### Fields

* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueknown::Bool`: Boolean flag for if the true solution is given.
* `u_analytic::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `timeseries`::AbstractArrayOrVoid`: u over time. Only saved if `save_timeseries=true`
  is specified in the solver.
* `t::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
  is specified in the solver.
* `timeseries_analytic`: If `save_timeseries=true`, saves the solution at each timestep.
* `prob::DEProblem`: Holds the problem object used to define the problem.
* `save_timeseries::Bool`: True if solver saved the extra timepoints.
* `appxtrue::Bool`: Boolean flag for if u_analytic was an approximation.

"""
type ODESolution <: AbstractODESolution
  u#::AbstractArrayOrNumber
  trueknown::Bool
  u_analytic#::AbstractArrayOrNumber
  errors#::Dict{}
  timeseries::AbstractArrayOrVoid
  t::AbstractArrayOrVoid
  timeseries_analytic::AbstractArrayOrVoid
  appxtrue::Bool
  save_timeseries::Bool
  k#::uType
  prob#
  alg
  interp::Function
  dense::Bool
  function ODESolution(u,prob,alg;timeseries=[],timeseries_analytic=[],t=[],k=[],saveat=[])
    save_timeseries = length(timeseries) > 2
    trueknown = false
    dense = k != []
    saveat_idxs = find((x)->x∈saveat,t)
    non_saveat_idxs = symdiff(1:length(t),saveat_idxs)

    t_nosaveat = @view t[non_saveat_idxs]
    timeseries_nosaveat = @view timeseries[non_saveat_idxs]
    if dense # dense
      if !prob.isinplace && typeof(u)<:AbstractArray
        f! = (t,u,du) -> (du[:] = prob.f(t,u))
      else
        f! = prob.f
      end
      interp = (tvals) -> ode_interpolation(tvals,t_nosaveat,timeseries_nosaveat,k,alg,f!)
    else
      interp = (tvals) -> nothing
    end
    return(new(u,trueknown,nothing,Dict(),timeseries,t,timeseries_analytic,false,save_timeseries,k,prob,alg,interp,dense))
  end
  function ODESolution(u,u_analytic,prob,alg;timeseries=[],timeseries_analytic=[],
           t=[],k=[],saveat=[],timeseries_errors=true,dense_errors=true)
    save_timeseries = length(timeseries) > 2
    trueknown = true

    dense = length(k)>1
    saveat_idxs = find((x)->x∈saveat,t)
    t_nosaveat = view(t,symdiff(1:length(t),saveat_idxs))
    timeseries_nosaveat = view(timeseries,symdiff(1:length(t),saveat_idxs))
    if dense # dense
      if !prob.isinplace && typeof(u)<:AbstractArray
        f! = (t,u,du) -> (du[:] = prob.f(t,u))
      else
        f! = prob.f
      end
      interp = (tvals) -> ode_interpolation(tvals,t_nosaveat,timeseries_nosaveat,k,alg,f!)
    else
      interp = (tvals) -> nothing
    end

    errors = Dict{Symbol,eltype(u)}()
    errors[:final] = mean(abs.(u-u_analytic))

    if save_timeseries && timeseries_errors
      errors[:l∞] = maximum(vecvecapply((x)->abs.(x),timeseries-timeseries_analytic))
      errors[:l2] = sqrt(mean(vecvecapply((x)->float(x).^2,timeseries-timeseries_analytic)))
      if dense && dense_errors
        densetimes = collect(linspace(t[1],t[end],100))
        interp_u = interp(densetimes)
        interp_analytic = [prob.analytic(t,timeseries[1]) for t in densetimes]
        errors[:L∞] = maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic))
        errors[:L2] = sqrt(mean(vecvecapply((x)->float(x).^2,interp_u-interp_analytic)))
      end
    end
    return(new(u,trueknown,u_analytic,errors,timeseries,t,timeseries_analytic,false,save_timeseries,k,prob,alg,interp,dense))
  end
end

(sol::ODESolution)(t) = sol.interp(t)
