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

# Despecializing levels erase callbacks already on the problem; skip when absent
# so merge_problem_kwargs does not hit the type-unstable CallbackSet merge path.
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
@test despecialized_solved === despecialized_default_problem
@test !haskey(despecialized_solved.kwargs, :callback)

if VERSION >= v"1.12"
    despecialized_cb = ContinuousCallback(
        (u, t, integrator) -> u - 0.5, integrator -> nothing
    )
    despecialized_with_cb = ODEProblem{false, SciMLBase.AutoSpecialize}(
        despecialized_default_f, 1.0, (0.0, 1.0), 2.0; callback = despecialized_cb
    )
    despecialized_cb_solved, = solve(despecialized_with_cb; wrap = Val(false))
    @test despecialized_cb_solved.kwargs[:callback] isa
        SciMLBase.CallbackSet{Vector{Any}, Vector{Any}}
    @test only(despecialized_cb_solved.kwargs[:callback].continuous_callbacks) ===
        despecialized_cb

    f_alloc!(du, u, p, t) = (du[1] = -u[1]; nothing)
    prob_alloc = ODEProblem{true, SciMLBase.AutoSpecialize}(f_alloc!, [1.0], (0.0, 1.0))
    cb_alloc = ContinuousCallback((u, t, i) -> u[1], i -> nothing)
    dcb_alloc = DiscreteCallback((u, t, i) -> false, i -> nothing)
    cbs_alloc = CallbackSet(cb_alloc, dcb_alloc)
    concrete_alloc = DiffEqBase.get_concrete_problem(prob_alloc, true)
    @test !haskey(concrete_alloc.kwargs, :callback)
    merge_once() = DiffEqBase.merge_problem_kwargs(concrete_alloc; callback = cbs_alloc)
    merge_once()
    @test (@allocated merge_once()) < 1500
end

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
