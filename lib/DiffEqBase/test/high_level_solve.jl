using DiffEqBase, Test
using Distributions
import SciMLBase

@test DiffEqBase.promote_tspan((0.0, 1.0)) == (0.0, 1.0)
@test DiffEqBase.promote_tspan((0, 1.0)) == (0.0, 1.0)
@test DiffEqBase.promote_tspan(1.0) == (0.0, 1.0)
@test DiffEqBase.promote_tspan(nothing) == (nothing, nothing)
@test DiffEqBase.promote_tspan(Real[0, 1.0]) == (0.0, 1.0)

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/1776
# promote_tspan(u0, p, tspan, prob, kwargs)
@test DiffEqBase.promote_tspan((0, 1)) == (0, 1)
@test DiffEqBase.promote_tspan(nothing, nothing, (0, 1), nothing, (dt = 1,)) == (0, 1)
@test DiffEqBase.promote_tspan(nothing, nothing, (0, 1), nothing, (dt = 1 / 2,)) ==
    (0.0, 1.0)

prob = ODEProblem((u, p, t) -> u, (p, t0) -> p[1], (p) -> (0.0, p[2]), (2.0, 1.0))
prob2 = DiffEqBase.get_concrete_problem(prob, true)

@test prob2.u0 == 2.0
@test prob2.tspan == (0.0, 1.0)

prob = ODEProblem((u, p, t) -> u, (p, t) -> Normal(p, 1), (0.0, 1.0), 1.0)
prob2 = DiffEqBase.get_concrete_problem(prob, true)
@test typeof(prob2.u0) == Float64

kwargs(; kw...) = kw
prob = ODEProblem((u, p, t) -> u, 1.0, nothing)
prob2 = DiffEqBase.get_concrete_problem(prob, true, tspan = (1.2, 3.4))
@test prob2.tspan === (1.2, 3.4)

prob = ODEProblem((u, p, t) -> u, nothing, nothing)
prob2 = DiffEqBase.get_concrete_problem(prob, true, u0 = 1.01, tspan = (1.2, 3.4))
@test prob2.u0 === 1.01

prob = ODEProblem((u, p, t) -> u, 1.0, (0, 1))
prob2 = DiffEqBase.get_concrete_problem(prob, true)
@test prob2.tspan == (0.0, 1.0)

prob = DDEProblem(
    (u, h, p, t) -> -h(p, t - p[1]), (p, t0) -> p[2], (p, t) -> 0,
    (p) -> (0.0, p[3]), (1.0, 2.0, 3.0); constant_lags = (p) -> [p[1]]
)
prob2 = DiffEqBase.get_concrete_problem(prob, true)

@test prob2.u0 == 2.0
@test prob2.tspan == (0.0, 3.0)
@test prob2.constant_lags == [1.0]

struct PreparedDefaultAlgorithm <: SciMLBase.AbstractODEAlgorithm end

prepared_default_problem = ODEProblem{false, SciMLBase.FullSpecialize}(
    (u, p, t) -> p * u, 1.0, (0.0, 1.0), 2.0
)
const PreparedDefaultProblem = typeof(prepared_default_problem)

DiffEqBase.prepare_alg(::Nothing, u0, p, ::PreparedDefaultProblem) =
    PreparedDefaultAlgorithm()
SciMLBase.__solve(prob::PreparedDefaultProblem, ::PreparedDefaultAlgorithm; kwargs...) =
    (prob, :solve)
SciMLBase.__init(prob::PreparedDefaultProblem, ::PreparedDefaultAlgorithm; kwargs...) =
    (prob, :init)

solved_problem, solve_stage = solve(prepared_default_problem; wrap = Val(false))
initialized_problem, init_stage = init(prepared_default_problem)
@test solved_problem === prepared_default_problem
@test initialized_problem === prepared_default_problem
@test solve_stage === :solve
@test init_stage === :init

# The despecializing levels store callbacks in type-erased vectors, so
# `get_concrete_problem` gives every such problem a `callback` kwarg -- including one
# built without any callback -- to keep the solver type constant as callbacks change.
despecialized_default_f(u, p, t) = p * u
despecialized_default_problem = ODEProblem{false, SciMLBase.AutoSpecialize}(
    despecialized_default_f, 1.0, (0.0, 1.0), 2.0
)
const DespecializedDefaultProblem = SciMLBase.ODEProblem{
    <:Any, <:Any, <:Any, <:Any,
    <:SciMLBase.ODEFunction{<:Any, <:Any, typeof(despecialized_default_f)},
}

DiffEqBase.prepare_alg(::Nothing, u0, p, ::DespecializedDefaultProblem) =
    PreparedDefaultAlgorithm()
SciMLBase.__solve(
    prob::DespecializedDefaultProblem, ::PreparedDefaultAlgorithm; kwargs...
) = (prob, :solve)

despecialized_solved, despecialized_stage = solve(
    despecialized_default_problem; wrap = Val(false)
)
@test despecialized_stage === :solve
@test despecialized_solved !== despecialized_default_problem
@test despecialized_solved.f === despecialized_default_problem.f
@test despecialized_solved.u0 == despecialized_default_problem.u0
@test despecialized_solved.tspan == despecialized_default_problem.tspan
@test despecialized_solved.kwargs[:callback] isa
    SciMLBase.CallbackSet{Vector{Any}, Vector{Any}}
@test isempty(despecialized_solved.kwargs[:callback].continuous_callbacks)
@test isempty(despecialized_solved.kwargs[:callback].discrete_callbacks)

# Problems that `ConstructionBase.setproperties` cannot rebuild are left untouched.
rode_problem = RODEProblem((u, p, t, W) -> u + W, 1.0, (0.0, 1.0))
bv_problem = BVProblem(
    (du, u, p, t) -> (du[1] = u[2]; du[2] = -u[1]),
    (res, u, p, t) -> (res[1] = u[1][1]; res[2] = u[end][1] - 1),
    [0.0, 0.0], (0.0, 1.0)
)
for prob in (rode_problem, bv_problem)
    @test DiffEqBase.get_concrete_problem(prob, true) === prob
    @test !haskey(prob.kwargs, :callback)
end
