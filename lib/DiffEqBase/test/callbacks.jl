using DiffEqBase, Test

condition = function (u, t, integrator) # Event when event_f(u,t,k) == 0
    return t - 2.95
end

affect! = function (integrator)
    return integrator.u = integrator.u + 2
end

rootfind = true
save_positions = (true, true)
callback = ContinuousCallback(condition, affect!; save_positions = save_positions)

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
saving_callback = DiscreteCallback(condition, affect!; save_positions = save_positions)

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
find_first_integrator = EmptyIntegrator([1.0, 2.0], 1, 0.0)
vector_affect! = function (integrator, idx)
    return integrator.u = integrator.u + idx
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
