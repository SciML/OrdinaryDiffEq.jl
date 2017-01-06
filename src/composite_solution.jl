### Concrete Types

type ODECompositeSolution{uType,tType,rateType,P,A} <: AbstractODESolution
  u::uType
  t::tType
  k::rateType
  prob::P
  alg::A
  interp::Function
  alg_choice::Vector{Int}
  dense::Bool
  tslocation::Int
end
(sol::ODECompositeSolution)(t) = sol.interp(t)

type ODECompositeTestSolution{uType,uType2,uEltype,tType,rateType,P,A} <: AbstractODETestSolution
  u::uType
  u_analytic::uType2
  errors::Dict{Symbol,uEltype}
  t::tType
  k::rateType
  prob::P
  alg::A
  interp::Function
  alg_choice::Vector{Int}
  dense::Bool
  tslocation::Int
end
(sol::ODECompositeTestSolution)(t) = sol.interp(t)

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
