using DiffEqBase, Test

import ReverseDiff, SciMLBase

@testset "ReverseDiff callback root derivative" begin
    for (left, right) in ((0.0, 1.0), (1.0, 0.0), (0.5, 1.0), (1.0, 0.5), (0.0, 0.5), (0.5, 0.0)),
            side in (SciMLBase.LeftRootFind, SciMLBase.RightRootFind)
        x = [0.5, left, right]
        function event_time(p)
            q, lo, hi = p
            left_residual, right_residual = lo - q, hi - q
            function condition(t, _)
                t == lo && return left_residual
                t == hi && return right_residual
                return t - q
            end
            return DiffEqBase.find_root(condition, (lo, hi), side)
        end
        @test event_time(x) == x[1]
        @test ReverseDiff.gradient(event_time, x) ≈ [1.0, 0.0, 0.0] atol = 1.0e-6
    end
end

@testset "ReverseDiff callback slope and replay" begin
    function event_time(f, x, side)
        condition(t, _) = f(t, x)
        return DiffEqBase.find_root(condition, (x[2], x[3]), side)
    end

    sharp(t, x) = atan((t - x[1]) / 1.0e-6)
    function endpoint_residual(t, x)
        delta = t - x[1]
        return t == x[2] || t == x[3] ? 2delta : delta
    end
    for (condition, q, lo, hi) in (
                (sharp, 0.5, 0.49998, 0.50002),
                (endpoint_residual, 0.625, prevfloat(0.625, 4), nextfloat(0.625, 4)),
            ),
            (a, b) in ((lo, hi), (hi, lo)),
            side in (SciMLBase.LeftRootFind, SciMLBase.RightRootFind)
        x = [q, a, b]
        root(x) = event_time(condition, x, side)
        @test root(x) == q
        @test ReverseDiff.gradient(root, x) ≈ [1.0, 0.0, 0.0] atol = 1.0e-6
    end

    affine(t, x) = t - x[1]
    scaled(t, x) = (1 + x[1]) * (t - x[1])
    squared_root(x) = event_time(affine, x, SciMLBase.LeftRootFind)^2
    scaled_root(x) = event_time(scaled, x, SciMLBase.LeftRootFind)
    recording = [0.25, 0.0, 1.0]
    replay = [0.375, 0.0, 1.0]
    for (f, expected) in ((squared_root, [0.75, 0.0, 0.0]), (scaled_root, [1.0, 0.0, 0.0]))
        tape = ReverseDiff.GradientTape(f, recording)
        for trace in (tape, ReverseDiff.compile(tape))
            @test ReverseDiff.gradient!(trace, replay) ≈ expected atol = 1.0e-6
        end
    end
end

@testset "Live event derivative decisions" begin
    function squared_crossing(condition, x, halfwidth, side, reversed)
        q = x[1]
        w = halfwidth(x)
        lo, hi = q - w, q + w
        f(t, _) = condition(t - q, x)
        bracket = reversed ? (hi, lo) : (lo, hi)
        return DiffEqBase.find_root(f, bracket, side)^2
    end

    cubic(s, x) = s * (1 + x[2] * s^2)
    narrow(s, x) = s * (1 - log1p(-(s / (2x[2]))^2))

    for (f, halfwidth, recording, replay, expected) in (
                (cubic, x -> 0.25, [0.5, 1.0], [0.625, 1.0e14], [1.25, 0.0]),
                (narrow, last, [0.5, 0.25], [0.5, 2.0^-21], [1.0, 0.0]),
            ),
            side in (SciMLBase.LeftRootFind, SciMLBase.RightRootFind), reversed in (false, true)
        g(x) = squared_crossing(f, x, halfwidth, side, reversed)
        @test g(replay) == replay[1]^2
        @test ReverseDiff.gradient(g, replay) ≈ expected atol = 1.0e-6
        tape = ReverseDiff.GradientTape(g, recording)
        for trace in (tape, ReverseDiff.compile(tape))
            @test ReverseDiff.gradient!(trace, replay) ≈ expected atol = 1.0e-6
        end
    end
end

@testset "ReverseDiff roots in adjacent brackets" begin
    q = 0.5
    for gap in (1, 2), edge in (prevfloat(q, gap), nextfloat(q, gap)),
            bracket in ((q, edge), (edge, q)),
            side in (SciMLBase.LeftRootFind, SciMLBase.RightRootFind)
        x = [q, bracket...]
        function root(x)
            return DiffEqBase.find_root((t, _) -> t - x[1], (x[2], x[3]), side)
        end
        @test root(x) == q
        @test ReverseDiff.gradient(root, x) ≈ [1.0, 0.0, 0.0] atol = 1.0e-6
    end
end

struct GenericEventArray{T} <: AbstractVector{T}
    data::Vector{T}
end
Base.size(array::GenericEventArray) = size(array.data)
Base.getindex(array::GenericEventArray, index::Int) = array.data[index]
Base.setindex!(array::GenericEventArray, value, index::Int) = array.data[index] = value

@testset "generic vector callback event storage" begin
    @test DiffEqBase.is_event_occurrence(-1.0, 1.0)
    @test DiffEqBase.is_event_occurrence(1.0, -1.0)
    @test DiffEqBase.is_event_occurrence(1.0, 0.0)
    @test !DiffEqBase.is_event_occurrence(0.0, 1.0)
    @test !DiffEqBase.is_event_occurrence(1.0, 1.0)
    @test !DiffEqBase.is_event_occurrence(-1.0, -1.0)

    next_generic = GenericEventArray([1.0, -1.0, 1.0, -1.0])
    prev_generic = GenericEventArray([-1.0, 1.0, 0.0, -1.0])
    @test DiffEqBase.findall_events!(next_generic, prev_generic)
    @test next_generic == [1.0, 1.0, 0.0, 0.0]

    next_array = [1.0, -1.0, 1.0, -1.0]
    prev_array = [-1.0, 1.0, 0.0, -1.0]
    @test DiffEqBase.findall_events!(next_array, prev_array)
    @test next_array == [1.0, 1.0, 0.0, 0.0]

    next_parent = [1.0, -1.0, 1.0, -1.0, 5.0]
    prev_parent = [-1.0, 1.0, 0.0, -1.0, 1.0]
    @test DiffEqBase.findall_events!(view(next_parent, 1:4), view(prev_parent, 1:4))
    @test next_parent == [1.0, 1.0, 0.0, 0.0, 5.0]
end


condition = function (u, t, integrator) # Event when event_f(u,t,k) == 0
    return t - 2.95
end

affect! = function (integrator)
    return integrator.u = integrator.u + 2
end

rootfind = true
save_positions = (true, true)
callback = ContinuousCallback(condition, affect!; save_positions)

cbs = CallbackSet(nothing)
@test typeof(cbs.discrete_callbacks) <: Tuple
@test typeof(cbs.continuous_callbacks) <: Tuple
cbs = CallbackSet(callback, nothing)
@test typeof(cbs.discrete_callbacks) <: Tuple
@test typeof(cbs.continuous_callbacks) <: Tuple
cbs = CallbackSet(callback, CallbackSet())
@test typeof(cbs.discrete_callbacks) <: Tuple
@test typeof(cbs.continuous_callbacks) <: Tuple

condition = function (integrator)
    return true
end
affect! = function (integrator) end
save_positions = (true, false)
saving_callback = DiscreteCallback(condition, affect!; save_positions)

cbs1 = CallbackSet(callback, saving_callback)

@test length(cbs1.discrete_callbacks) == 1
@test length(cbs1.continuous_callbacks) == 1

cbs2 = CallbackSet(callback)
@test length(cbs2.continuous_callbacks) == 1
@test length(cbs2.discrete_callbacks) == 0

cbs3 = CallbackSet(saving_callback)
@test length(cbs3.discrete_callbacks) == 1
@test length(cbs3.continuous_callbacks) == 0

cbs4 = CallbackSet()
@test length(cbs4.discrete_callbacks) == 0
@test length(cbs4.continuous_callbacks) == 0

cbs5 = CallbackSet(cbs1, cbs2)

@test length(cbs5.discrete_callbacks) == 1
@test length(cbs5.continuous_callbacks) == 2

# For the purposes of this test, create a empty integrator type and
# override find_callback_time, since we don't actually care about testing
# the find callback time aspect, just the inference failure
mutable struct EmptyIntegrator
    u::Vector{Float64}
    tdir::Int
    last_event_error::Float64
    callback_cache::DiffEqBase.CallbackCache{Vector{Float64}, Vector{Float64}}
    event_last_time::Int
end
function DiffEqBase.find_callback_time(
        integrator::EmptyIntegrator,
        callback::ContinuousCallback, counter
    )
    return 1.0 + counter, 0.9 + counter, true, counter, 0.0
end
function DiffEqBase.find_callback_time(
        integrator::EmptyIntegrator,
        callback::VectorContinuousCallback, counter
    )
    return 1.0 + counter, 0.9 + counter, true, counter, 0.0
end
find_first_integrator = EmptyIntegrator(
    [1.0, 2.0], 1, 0.0,
    DiffEqBase.CallbackCache(2, Float64, Float64),
    0,
)
vector_affect! = function (integrator, events)
    for (idx, dir) in enumerate(events)
        iszero(dir) && continue
        integrator.u = integrator.u .+ idx
    end
    return
end

cond_1(u, t, integrator) = t - 1.0
cond_2(u, t, integrator) = t - 1.1
cond_3(u, t, integrator) = t - 1.2
cond_4(u, t, integrator) = t - 1.3
cond_5(u, t, integrator) = t - 1.4
cond_6(u, t, integrator) = t - 1.5
cond_7(u, t, integrator) = t - 1.6
cond_8(u, t, integrator) = t - 1.7
cond_9(u, t, integrator) = t - 1.8
cond_10(u, t, integrator) = t - 1.9
# Setup a lot of callbacks so the recursive inference failure happens
callbacks = (
    ContinuousCallback(cond_1, affect!),
    ContinuousCallback(cond_2, affect!),
    ContinuousCallback(cond_3, affect!),
    ContinuousCallback(cond_4, affect!),
    ContinuousCallback(cond_5, affect!),
    ContinuousCallback(cond_6, affect!),
    ContinuousCallback(cond_7, affect!),
    ContinuousCallback(cond_8, affect!),
    ContinuousCallback(cond_9, affect!),
    ContinuousCallback(cond_10, affect!),
    VectorContinuousCallback(cond_1, vector_affect!, 2),
    VectorContinuousCallback(cond_2, vector_affect!, 2),
    VectorContinuousCallback(cond_3, vector_affect!, 2),
    VectorContinuousCallback(cond_4, vector_affect!, 2),
    VectorContinuousCallback(cond_5, vector_affect!, 2),
    VectorContinuousCallback(cond_6, vector_affect!, 2),
);
function test_find_first_callback(callbacks, int)
    return @timed(DiffEqBase.find_first_continuous_callback(int, callbacks...))
end
test_find_first_callback(callbacks, find_first_integrator);
@test test_find_first_callback(callbacks, find_first_integrator).bytes == 0

# https://github.com/SciML/DiffEqBase.jl/issues/1233
@testset "Inexact rootfinding" begin
    # function with irrational root (sqrt(2))
    irrational_f(x, p = 0.0) = x^2 - 2

    # Forward integration
    is_forward = true
    tspan = (1.0, 2.0)
    before = DiffEqBase.find_root(irrational_f, tspan, SciMLBase.LeftRootFind)
    after = DiffEqBase.find_root(irrational_f, tspan, SciMLBase.RightRootFind)
    @test irrational_f(before) < 0.0
    @test irrational_f(after) > 0.0
    @test nextfloat(before) == after

    # Backward integration
    is_forward = false
    tspan = (2.0, 1.0)
    before = DiffEqBase.find_root(irrational_f, tspan, SciMLBase.LeftRootFind)
    after = DiffEqBase.find_root(irrational_f, tspan, SciMLBase.RightRootFind)
    @test irrational_f(before) > 0.0
    @test irrational_f(after) < 0.0
    @test nextfloat(after) == before
end
