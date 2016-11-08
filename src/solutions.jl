"""
`ODESolution`

Holds the data for the solution to an ODE problem.

### Fields

* `u::Array{Float64}`: The solution
* `u_analytic::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `t::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
  is specified in the solver.
* `prob::DEProblem`: Holds the problem object used to define the problem.

"""
type ODESolution <: AbstractODESolution
  u#::AbstractArrayOrNumber
  u_analytic#::AbstractArrayOrNumber
  errors#::Dict{}
  t::AbstractArrayOrVoid
  k#::uType
  prob#
  alg
  interp::Function
  dense::Bool
end

function ODESolution{uType,tType,isinplace}(t,u,
        prob::AbstractODEProblem{uType,tType,Val{isinplace}},
        alg;u_analytic=[],
        k=[],saveat=[],timeseries_errors=true,dense_errors=true)
  save_timeseries = length(u) > 2

  dense = length(k)>1
  saveat_idxs = find((x)->x∈saveat,t)
  t_nosaveat = view(t,symdiff(1:length(t),saveat_idxs))
  u_nosaveat = view(u,symdiff(1:length(t),saveat_idxs))
  if dense # dense
    if !isinplace && typeof(u[1])<:AbstractArray
      f! = (t,u,du) -> (du[:] = prob.f(t,u))
    else
      f! = prob.f
    end
    interp = (tvals) -> ode_interpolation(tvals,t_nosaveat,u_nosaveat,k,alg,f!)
  else
    interp = (tvals) -> nothing
  end

  errors = Dict{Symbol,eltype(u[1])}()
  if !isempty(u_analytic)
    errors[:final] = mean(abs.(u[end]-u_analytic[end]))

    if save_timeseries && timeseries_errors
      errors[:l∞] = maximum(vecvecapply((x)->abs.(x),u-u_analytic))
      errors[:l2] = sqrt(mean(vecvecapply((x)->float(x).^2,u-u_analytic)))
      if dense && dense_errors
        densetimes = collect(linspace(t[1],t[end],100))
        interp_u = interp(densetimes)
        interp_analytic = [prob.analytic(t,u[1]) for t in densetimes]
        errors[:L∞] = maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic))
        errors[:L2] = sqrt(mean(vecvecapply((x)->float(x).^2,interp_u-interp_analytic)))
      end
    end
  end
  return(ODESolution(u,u_analytic,errors,t,k,prob,alg,interp,dense))
end

(sol::ODESolution)(t) = sol.interp(t)
