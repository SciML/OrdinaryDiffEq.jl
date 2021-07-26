struct ODECompositeSolution{T,N,uType,uType2,EType,tType,rateType,P,A,IType,DE} <: DiffEqBase.AbstractODESolution{T,N,uType}
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
  destats::DE
  retcode::Symbol
  residual::Function
end
(sol::ODECompositeSolution)(t,deriv::Type=Val{0};idxs=nothing,continuity=:left) = sol.interp(t,idxs,deriv,sol.prob.p,continuity)
(sol::ODECompositeSolution)(v,t,deriv::Type=Val{0};idxs=nothing,continuity=:left) = sol.interp(v,t,idxs,deriv,sol.prob.p,continuity)

function DiffEqBase.build_solution(
        prob::Union{DiffEqBase.AbstractODEProblem,DiffEqBase.AbstractDDEProblem},
        alg::OrdinaryDiffEqCompositeAlgorithm,t,u;
        timeseries_errors=length(u)>2,
        dense=false,dense_errors=dense,
        calculate_error = true,
        k=[],
        du=[],
        interp = !isempty(du) ? 
          DiffEqBase.HermiteInterpolation(t,u,du) : DiffEqBase.LinearInterpolation(t,u),
        alg_choice=[1], retcode = :Default, destats=DiffEqBase.DEStats(), kwargs...)

  T = eltype(eltype(u))
  if typeof(prob.u0) <: Tuple
    N = length((size(ArrayPartition(prob.u0))..., length(u)))
  else
    N = length((size(prob.u0)..., length(u)))
  end

  f = prob.f

  residual(t) = interp(t,Val{1}) - f(interp(t), prob.p, t)
  if DiffEqBase.has_analytic(f)
    if !(typeof(prob.u0) <: Tuple)
      u_analytic = Vector{typeof(prob.u0)}(undef, 0)
      errors = Dict{Symbol,real(eltype(prob.u0))}()
    else
      u_analytic = Vector{typeof(ArrayPartition(prob.u0))}(undef, 0)
      errors = Dict{Symbol,real(eltype(prob.u0[1]))}()
    end

    sol = ODECompositeSolution{T,N,typeof(u),typeof(u_analytic),typeof(errors),typeof(t),typeof(k),
                               typeof(prob),typeof(alg),typeof(interp),typeof(destats)}(u,u_analytic,errors,t,k,prob,alg,interp,alg_choice,dense,0,destats,retcode,residual)
    if calculate_error
      DiffEqBase.calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
    end
    return sol
  else
    return ODECompositeSolution{T,N,typeof(u),typeof(nothing),typeof(nothing),typeof(t),typeof(k),
                                typeof(prob),typeof(alg),typeof(interp),typeof(destats)}(u,nothing,nothing,t,k,prob,alg,interp,alg_choice,dense,0,destats,retcode,residual)
  end
end

function DiffEqBase.solution_new_retcode(sol::ODECompositeSolution{T,N,uType,uType2,EType,tType,rateType,P,A,IType,DE},retcode) where {T,N,uType,uType2,EType,tType,rateType,P,A,IType,DE}
  residual(t) = sol(t,Val{1}) - sol.prob.f(sol(t), sol.prob.p, t)
  ODECompositeSolution{T,N,uType,uType2,EType,tType,rateType,P,A,IType,DE}(
                     sol.u,sol.u_analytic,sol.errors,sol.t,sol.k,sol.prob,
                     sol.alg,sol.interp,sol.alg_choice,sol.dense,sol.tslocation,
                     sol.destats,retcode,residual)
end

function DiffEqBase.solution_new_tslocation(sol::ODECompositeSolution{
   T,N,uType,uType2,EType,tType,rateType,P,A,IType,DE},tslocation) where {T,N,uType,uType2,EType,tType,rateType,P,A,IType,DE}
   residual(t) = sol(t,Val{1}) - sol.prob.f(sol(t), sol.prob.p, t)
   ODECompositeSolution{T,N,uType,uType2,EType,tType,rateType,P,A,IType,DE}(
                      sol.u,sol.u_analytic,sol.errors,sol.t,sol.k,sol.prob,
                      sol.alg,sol.interp,sol.alg_choice,sol.dense,tslocation,
                      sol.destats,sol.retcode,residual)
end

function DiffEqBase.sensitivity_solution(sol::ODECompositeSolution,u,t)
    T = eltype(eltype(u))
    N = length((size(sol.prob.u0)..., length(u)))
    interp = if !sol.dense
      sol.interp
    else
      DiffEqBase.SensitivityInterpolation(t,u)
    end

    residual(t) = sol(t,Val{1}) - sol.prob.f(sol(t), sol.prob.p, t)
    ODECompositeSolution{T,N,typeof(u),typeof(sol.u_analytic),typeof(sol.errors),
                typeof(t),Nothing,typeof(sol.prob),typeof(sol.alg),
                typeof(interp),typeof(sol.destats)}(
                u,sol.u_analytic,sol.errors,t,nothing,sol.prob,
                sol.alg,interp,sol.alg_choice,
                sol.dense,sol.tslocation,
                sol.destats,sol.retcode,residual)
end
