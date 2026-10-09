using DiffEqBase, Test, ForwardDiff
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

@testset "Integer tspan promotion for tstops" begin
    span = (0, 1)
    promote_span(kw) = DiffEqBase.promote_tspan(nothing, nothing, span, nothing, kw)
    @test promote_span((tstops = [0.33, 0.66, 1.0],)) === (0.0, 1.0)
    @test promote_span((tstops = (0.33, 0.66, 1.0),)) === (0.0, 1.0)
    @test promote_span((tstops = 0.33,)) === (0.0, 1.0)
    @test promote_span((tstops = Any[0, 0.33, 1],)) === (0.0, 1.0)
    @test promote_span((tstops = (0, 0.33, 1),)) === (0.0, 1.0)
    @test promote_span((tstops = [1 // 3, 2 // 3],)) === (0 // 1, 1 // 1)
    @test promote_span((tstops = [0, 1],)) === span
    @test promote_span((tstops = Float64[],)) === span
    @test promote_span((tstops = (),)) === span
    @test promote_span((tstops = (p, tspan) -> [0.33],)) === span
    @test promote_span((dt = 1 // 3, tstops = [1 // 3, 2 // 3])) === (0 // 1, 1 // 1)
    @test promote_span((dt = 0.1, tstops = [0.33])) === (0.0, 1.0)
    large_span = (2^53, 2^53 + 1)
    for stops in ([0.5], (Float64(2^53),), Any[Float64(2^53 + 2)])
        @test DiffEqBase.promote_tspan(nothing, nothing, large_span, nothing, (tstops = stops,)) === large_span
    end
    @test promote_span((tstops = [-0.5, 0.0, 1.0, 1.5],)) === span
    @test DiffEqBase.promote_tspan(nothing, nothing, (0, 2), nothing, (tstops = [1.0],)) === (0, 2)
    @test DiffEqBase.promote_tspan(nothing, nothing, (0, 2), nothing, (tstops = [1 // 1],)) === (0, 2)
    @test promote_span((tstops = [0.0, 0.5],)) === (0.0, 1.0)
    @test promote_span((tstops = (0.5f0, 0.25),)) === (0.0, 1.0)
    @test DiffEqBase.promote_tspan(
        nothing, nothing, (0.0f0, 1.0f0), nothing, (tstops = [0.33, 0.66],)
    ) === (0.0f0, 1.0f0)
end

@testset "Problem time options before promotion" begin
    tstops = [0.33, 0.66, 1.0]
    ode = ODEProblem((du, u, p, t) -> du .= u, [1.0], (0, 1); tstops)
    dae = DAEProblem((res, du, u, p, t) -> res .= du .- u, [1.0], [1.0], (0, 1); tstops)
    dde = DDEProblem((du, u, h, p, t) -> du .= u, [1.0], (p, t) -> [1.0], (0, 1); tstops)
    for prob in (ode, dae, dde)
        @test DiffEqBase.get_concrete_problem(prob, false).tspan === (0.0, 1.0)
        @test DiffEqBase.get_concrete_problem(prob, false; tstops = [1]).tspan === (0, 1)
        @test DiffEqBase.get_concrete_problem(prob, false; tstops = [1 // 2]).tspan === (0 // 1, 1 // 1)
        @test DiffEqBase.get_concrete_problem(prob, false; tstops = Float64[]).tspan === (0, 1)
        @test DiffEqBase.get_concrete_problem(prob, false; tstops = ()).tspan === (0, 1)
    end
    dt32 = remake(ode; dt = 0.1f0)
    @test DiffEqBase.get_concrete_problem(dt32, false).tspan === (0.0f0, 1.0f0)
    dt64 = remake(ode; dt = 0.1)
    @test DiffEqBase.get_concrete_problem(dt64, false; dt = 0.1f0).tspan === (0.0f0, 1.0f0)
    span32 = remake(ode; tspan = (0.0f0, 1.0f0))
    @test DiffEqBase.get_concrete_problem(span32, false).tspan === (0.0f0, 1.0f0)
    tag = ForwardDiff.Tag(identity, Float32)
    dual_u0 = [ForwardDiff.Dual{typeof(tag)}(1.0f0, 1.0f0)]
    dual_prob = remake(ode; u0 = dual_u0)
    @test DiffEqBase.get_concrete_problem(dual_prob, false).tspan === (0.0, 1.0)
    complex_prob = remake(dt64; u0 = complex.(dual_u0))
    complex_span = DiffEqBase.get_concrete_problem(complex_prob, false).tspan
    @test typeof(ForwardDiff.value(complex_span[1])) === Float64
    complex_span32 = DiffEqBase.get_concrete_problem(complex_prob, false; dt = 0.1f0).tspan
    @test typeof(ForwardDiff.value(complex_span32[1])) === Float32
    dual_span = (ForwardDiff.Dual{typeof(tag)}(0.0f0, 0.0f0), ForwardDiff.Dual{typeof(tag)}(1.0f0, 0.0f0))
    dual_span_prob = remake(dt64; u0 = dual_u0, tspan = dual_span)
    @test typeof(ForwardDiff.value(DiffEqBase.get_concrete_problem(dual_span_prob, false).tspan[1])) === Float64
    @test typeof(ForwardDiff.value(DiffEqBase.get_concrete_problem(dual_span_prob, false; dt = 0.1f0).tspan[1])) === Float32
end
