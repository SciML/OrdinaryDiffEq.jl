### Concrete Types

type ODECompositeSolution{uType,tType,rateType,P,A,IType} <: AbstractODESolution
  u::uType
  t::tType
  k::rateType
  prob::P
  alg::A
  interp::IType
  alg_choice::Vector{Int}
  dense::Bool
  tslocation::Int
end
(sol::ODECompositeSolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::ODECompositeSolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

type ODECompositeTestSolution{uType,uType2,uEltype,tType,rateType,P,A,IType} <: AbstractODETestSolution
  u::uType
  u_analytic::uType2
  errors::Dict{Symbol,uEltype}
  t::tType
  k::rateType
  prob::P
  alg::A
  interp::IType
  alg_choice::Vector{Int}
  dense::Bool
  tslocation::Int
end
(sol::ODECompositeTestSolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::ODECompositeTestSolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

function build_solution{uType,tType,isinplace}(
        prob::AbstractODEProblem{uType,tType,isinplace},
        alg::OrdinaryDiffEqCompositeAlgorithm,t,u;dense=false,alg_choice=[1],
        k=[],interp = (tvals) -> nothing,kwargs...)
  ODECompositeSolution(u,t,k,prob,alg,interp,alg_choice,dense,0)
end

function build_solution{uType,tType,isinplace}(
        prob::AbstractODETestProblem{uType,tType,isinplace},
        alg::OrdinaryDiffEqCompositeAlgorithm,t,u,;dense=false,
        k=[],interp = (tvals) -> nothing,alg_choice=[1],
        timeseries_errors=true,dense_errors=true,
        calculate_error = true,kwargs...)
  u_analytic = Vector{uType}(0)
  errors = Dict{Symbol,eltype(u[1])}()
  sol = ODECompositeTestSolution(u,u_analytic,errors,t,k,prob,alg,interp,alg_choice,dense,0)
  if calculate_error
    calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
  end
  sol
end
