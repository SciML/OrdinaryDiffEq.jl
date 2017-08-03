mutable struct ODECompositeSolution{T,N,uType,uType2,EType,tType,rateType,P,A,IType} <: AbstractODESolution{T,N}
  u::uType
  u_analytic::uType2
  errors::EType
  t::tType
  k::rateType
  prob::P
  alg::A
  interp::IType
  alg_choice::Vector{Int}
  dense::Bool
  tslocation::Int
  retcode::Symbol
end
(sol::ODECompositeSolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::ODECompositeSolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

function build_solution{uType,tType,isinplace}(
        prob::AbstractODEProblem{uType,tType,isinplace},
        alg::OrdinaryDiffEqCompositeAlgorithm,t,u;dense=false,alg_choice=[1],
        k=[],interp = (tvals) -> nothing,
        timeseries_errors=true,dense_errors=true,
        calculate_error = true, retcode = :Default, kwargs...)

  T = eltype(eltype(u))
  if typeof(prob.u0) <: Tuple
    N = length((size(ArrayPartition(prob.u0))..., length(u)))
  else
    N = length((size(prob.u0)..., length(u)))
  end

  if typeof(prob.f) <: Tuple
    f = prob.f[1]
  else
    f = prob.f
  end

  if has_analytic(f)
    u_analytic = Vector{uType}(0)
    errors = Dict{Symbol,eltype(u[1])}()
    sol = ODECompositeSolution{T,N,typeof(u),typeof(u_analytic),typeof(errors),typeof(t),typeof(k),
                       typeof(prob),typeof(alg),typeof(interp)}(u,u_analytic,errors,t,k,prob,alg,interp,alg_choice,dense,0,retcode)
    if calculate_error
      calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
    end
    return sol
  else
    return ODECompositeSolution{T,N,typeof(u),typeof(u_analytic),typeof(errors),typeof(t),typeof(k),
                       typeof(prob),typeof(alg),typeof(interp)}(u,nothing,nothing,t,k,prob,alg,interp,alg_choice,dense,0,retcode)
  end
end
