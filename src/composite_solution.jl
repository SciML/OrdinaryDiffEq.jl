struct ODECompositeSolution{T,N,uType,uType2,EType,tType,rateType,P,A,IType} <: DiffEqBase.AbstractODESolution{T,N}
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
(sol::ODECompositeSolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv,sol.prob.p)
(sol::ODECompositeSolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv,sol.prob.p)

function DiffEqBase.build_solution(
        prob::Union{DiffEqBase.AbstractODEProblem,DiffEqBase.AbstractDDEProblem},
        alg::OrdinaryDiffEqCompositeAlgorithm,t,u;
        timeseries_errors=length(u)>2,
        dense=false,dense_errors=dense,
        calculate_error = true,
        k=[],
        du=[],
        interp = !isempty(du) ? HermiteInterpolation(t,u,du) : LinearInterpolation(t,u),
        alg_choice=[1], retcode = :Default, kwargs...)

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

  if DiffEqBase.has_analytic(f)
    if !(typeof(prob.u0) <: Tuple)
      u_analytic = Vector{typeof(prob.u0)}(undef, 0)
      errors = Dict{Symbol,real(eltype(prob.u0))}()
    else
      u_analytic = Vector{typeof(ArrayPartition(prob.u0))}(undef, 0)
      errors = Dict{Symbol,real(eltype(prob.u0[1]))}()
    end

    sol = ODECompositeSolution{T,N,typeof(u),typeof(u_analytic),typeof(errors),typeof(t),typeof(k),
                       typeof(prob),typeof(alg),typeof(interp)}(u,u_analytic,errors,t,k,prob,alg,interp,alg_choice,dense,0,retcode)
    if calculate_error
      DiffEqBase.calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
    end
    return sol
  else
    return ODECompositeSolution{T,N,typeof(u),typeof(nothing),typeof(nothing),typeof(t),typeof(k),
                       typeof(prob),typeof(alg),typeof(interp)}(u,nothing,nothing,t,k,prob,alg,interp,alg_choice,dense,0,retcode)
  end
end

function DiffEqBase.solution_new_retcode(sol::ODECompositeSolution{T,N,uType,uType2,EType,tType,rateType,P,A,IType},retcode) where {T,N,uType,uType2,EType,tType,rateType,P,A,IType}
  ODECompositeSolution{T,N,uType,uType2,EType,tType,rateType,P,A,IType}(
                     sol.u,sol.u_analytic,sol.errors,sol.t,sol.k,sol.prob,
                     sol.alg,sol.interp,sol.alg_choice,sol.dense,sol.tslocation,
                     retcode)
 end

 function DiffEqBase.solution_new_tslocation(sol::ODECompositeSolution{
   T,N,uType,uType2,EType,tType,rateType,P,A,IType},tslocation) where {T,N,uType,uType2,EType,tType,rateType,P,A,IType}
   ODECompositeSolution{T,N,uType,uType2,EType,tType,rateType,P,A,IType}(
                      sol.u,sol.u_analytic,sol.errors,sol.t,sol.k,sol.prob,
                      sol.alg,sol.interp,sol.alg_choice,sol.dense,tslocation,
                      sol.retcode)
  end
