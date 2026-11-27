using OrdinaryDiffEq, SciMLBase, Test

function callback_erasure_decay!(du, u, p, t)
    du[1] = -u[1]
    return nothing
end

@testset "callback erasure opt-out" begin
    continuous = ContinuousCallback(
        (u, t, integrator) -> sin(10π * t + 0.5), integrator -> nothing;
        save_positions = (false, false)
    )
    discrete = DiscreteCallback(
        (u, t, integrator) -> true, integrator -> nothing;
        save_positions = (false, false)
    )
    for specialize in (SciMLBase.AutoSpecialize, SciMLBase.AutoDespecialize, SciMLBase.NoSpecialize)
        prob = ODEProblem{true, specialize}(callback_erasure_decay!, [1.0], (0.0, 1.0); callback = continuous)
        integ = init(prob, Tsit5(); callback = discrete, erase_callbacks = false)
        @test integ.opts.callback.continuous_callbacks isa Tuple
        @test integ.opts.callback.discrete_callbacks isa Tuple
        @test SciMLBase.successful_retcode(solve!(integ))
        @test SciMLBase.successful_retcode(solve(prob, Tsit5(); callback = discrete, erase_callbacks = false))
    end
    prob = ODEProblem{true, SciMLBase.AutoSpecialize}(
        callback_erasure_decay!, [1.0], (0.0, 1.0); callback = continuous,
        erase_callbacks = false
    )
    @test init(prob, Tsit5(); callback = discrete).opts.callback.continuous_callbacks isa Tuple
    if VERSION >= v"1.12"
        @test init(prob, Tsit5(); erase_callbacks = true).opts.callback.continuous_callbacks isa Vector{Any}
    end
end

function callback_allocations(tend, callback)
    prob = ODEProblem{true, SciMLBase.AutoSpecialize}(
        callback_erasure_decay!, [1.0], (0.0, tend)
    )
    solve(prob, Tsit5(); dtmax = 0.01, save_everystep = false, callback, erase_callbacks = false)
    return @allocations solve(
        prob, Tsit5(); dtmax = 0.01, save_everystep = false, callback,
        erase_callbacks = false
    )
end

@testset "typed callback allocation scaling" begin
    callbacks = (
        ContinuousCallback(
            (u, t, integrator) -> sin(10π * t + 0.5), integrator -> nothing;
            save_positions = (false, false)
        ),
        DiscreteCallback(
            (u, t, integrator) -> true, integrator -> nothing;
            save_positions = (false, false)
        ),
    )
    for callback in callbacks
        short = callback_allocations(1.0, callback)
        long = callback_allocations(10.0, callback)
        @test long - short < 200
    end
end

@testset "callback erasure sentinel" begin
    cb = ContinuousCallback((u, t, integrator) -> t - 0.5, integrator -> nothing)
    prob = ODEProblem(callback_erasure_decay!, [1.0], (0.0, 1.0); callback = cb, erase_callbacks = false)
    @test init(prob, Tsit5(); erase_callbacks = nothing).opts.callback.continuous_callbacks isa Tuple
    @test SciMLBase.successful_retcode(solve(prob, Tsit5(); erase_callbacks = nothing))
    unset = ODEProblem(callback_erasure_decay!, [1.0], (0.0, 1.0); erase_callbacks = nothing)
    @test SciMLBase.successful_retcode(solve(unset, Tsit5()))
    @test_throws ArgumentError solve(prob, Tsit5(); erase_callbacks = 0)
end

import DiffEqBase

mutable struct DeferredCallbackWrapper <: SciMLBase.AbstractJumpProblem{Nothing, Nothing}
    prob::Any
    kwargs::Base.Pairs
    DeferredCallbackWrapper() = new()
end

function DeferredCallbackWrapper(kwargs)
    wrapper = DeferredCallbackWrapper()
    wrapper.kwargs = pairs(kwargs)
    return wrapper
end

SciMLBase.__init(::DeferredCallbackWrapper, ::Tsit5; save_everystep = true, kwargs...) = save_everystep
function SciMLBase.__solve(prob::DeferredCallbackWrapper, ::Tsit5; kwargs...)
    return DiffEqBase.merge_problem_kwargs(prob; kwargs...)[:save_everystep]
end

@testset "deferred wrapper keyword merging" begin
    prob = DeferredCallbackWrapper((;))
    for extra in ((;), (; erase_callbacks = nothing))
        @test !DiffEqBase.merge_problem_kwargs(prob; save_everystep = false, extra...)[:save_everystep]
        @test !init(prob, Tsit5(); save_everystep = false, extra...)
        @test !solve(prob, Tsit5(); save_everystep = false, extra...)
    end
end

mutable struct DeferredForwardingWrapper <: SciMLBase.AbstractJumpProblem{Nothing, Nothing}
    prob::Any
    kwargs::Base.Pairs
    default_erase::Bool
    function DeferredForwardingWrapper(kwargs, default_erase)
        wrapper = new()
        wrapper.kwargs = pairs(kwargs)
        wrapper.default_erase = default_erase
        return wrapper
    end
end

function SciMLBase.__init(wrapper::DeferredForwardingWrapper, alg::Tsit5; kwargs...)
    merged = DiffEqBase.merge_problem_kwargs(wrapper; kwargs...)
    wrapper.prob = ODEProblem{true, SciMLBase.AutoSpecialize}(
        callback_erasure_decay!, [1.0], (0.0, 1.0)
    )
    return init(wrapper.prob, alg; erase_callbacks = wrapper.default_erase, merged...)
end

function SciMLBase.__solve(wrapper::DeferredForwardingWrapper, alg::Tsit5; kwargs...)
    return SciMLBase.__init(wrapper, alg; kwargs...)
end

@testset "deferred wrapper forwards resolved erasure" begin
    callbacks = CallbackSet(
        ContinuousCallback((u, t, integrator) -> t - 0.5, integrator -> nothing),
        DiscreteCallback((u, t, integrator) -> false, integrator -> nothing)
    )
    cases = (
        ((;), false, false),
        ((; erase_callbacks = false), nothing, false),
        ((; erase_callbacks = false), true, true),
    )
    for (stored, setting, resolved) in cases
        wrapper = () -> DeferredForwardingWrapper(stored, !resolved)
        @test DiffEqBase.merge_problem_kwargs(wrapper(); erase_callbacks = setting)[:erase_callbacks] === resolved
        for run in (init, solve)
            cb = () -> run(wrapper(), Tsit5(); callback = callbacks, erase_callbacks = setting).opts.callback
            if resolved
                VERSION >= v"1.12" && @test cb().continuous_callbacks isa Vector{Any}
            else
                @test cb().continuous_callbacks isa Tuple
                @test cb().discrete_callbacks isa Tuple
            end
        end
    end
end

struct CallbackWrapper{P, C, K} <: SciMLBase.AbstractJumpProblem{P, Nothing}
    prob::P
    callback::C
    kwargs::K
    default_erase::Bool
end

function SciMLBase.__init(prob::CallbackWrapper, alg::Tsit5; kwargs...)
    merged = DiffEqBase.merge_problem_kwargs(prob; kwargs...)
    return init(prob.prob, alg; erase_callbacks = prob.default_erase, callback = prob.callback, merged...)
end

SciMLBase.__solve(prob::CallbackWrapper, alg::Tsit5; kwargs...) = SciMLBase.__init(prob, alg; kwargs...)

struct CallbackTestSenseAlg end
function SciMLBase._concrete_solve_adjoint(prob, alg, ::CallbackTestSenseAlg, u0, p, originator; kwargs...)
    return init(prob, alg; kwargs...)
end
function SciMLBase._concrete_solve_forward(prob, alg, ::CallbackTestSenseAlg, u0, p, originator; kwargs...)
    return init(prob, alg; kwargs...)
end

@testset "callback erasure forwarding" begin
    continuous = ContinuousCallback((u, t, integrator) -> t - 0.5, integrator -> nothing)
    discrete = DiscreteCallback((u, t, integrator) -> false, integrator -> nothing)
    callbacks = CallbackSet(continuous, discrete)
    prob = ODEProblem(callback_erasure_decay!, [1.0], (0.0, 1.0))
    wrapper = CallbackWrapper(prob, callbacks, pairs((;)), true)
    for run in (init, solve)
        cb = run(wrapper, Tsit5(); erase_callbacks = false).opts.callback
        @test cb.continuous_callbacks isa Tuple
        @test cb.discrete_callbacks isa Tuple
    end
    stored = CallbackWrapper(prob, callbacks, pairs((; erase_callbacks = false)), true)
    @test init(stored, Tsit5()).opts.callback.continuous_callbacks isa Tuple
    if VERSION >= v"1.12"
        @test init(CallbackWrapper(prob, callbacks, pairs((;)), false), Tsit5(); erase_callbacks = true).opts.callback.continuous_callbacks isa Vector{Any}
    end
    prob_cb = ODEProblem(callback_erasure_decay!, [1.0], (0.0, 1.0); callback = callbacks)
    prob_optout = ODEProblem(callback_erasure_decay!, [1.0], (0.0, 1.0); callback = callbacks, erase_callbacks = false)
    for run in (DiffEqBase._solve_adjoint, DiffEqBase._solve_forward)
        cb = run(prob_cb, CallbackTestSenseAlg(), prob_cb.u0, prob_cb.p, nothing, Tsit5(); erase_callbacks = false).opts.callback
        @test cb.continuous_callbacks isa Tuple
        @test cb.discrete_callbacks isa Tuple
        if VERSION >= v"1.12"
            cb = run(prob_optout, CallbackTestSenseAlg(), prob_optout.u0, prob_optout.p, nothing, Tsit5(); erase_callbacks = true).opts.callback
            @test cb.continuous_callbacks isa Vector{Any}
        end
    end
end

struct PropertyForwardingWrapper{K} <: SciMLBase.AbstractJumpProblem{Nothing, Nothing}
    prob::Nothing
    kwargs::K
end

function Base.getproperty(wrapper::PropertyForwardingWrapper, name::Symbol)
    name === :prob || return getfield(wrapper, name)
    return ODEProblem{true, SciMLBase.AutoSpecialize}(callback_erasure_decay!, [1.0], (0.0, 1.0))
end

function SciMLBase.__init(wrapper::PropertyForwardingWrapper, alg::Tsit5; kwargs...)
    return init(wrapper.prob, alg; DiffEqBase.merge_problem_kwargs(wrapper; kwargs...)...)
end

function SciMLBase.__solve(wrapper::PropertyForwardingWrapper, alg::Tsit5; kwargs...)
    return SciMLBase.__init(wrapper, alg; kwargs...)
end

@testset "property-backed wrapper forwards erasure" begin
    callbacks = CallbackSet(
        ContinuousCallback((u, t, integrator) -> t - 0.5, integrator -> (integrator.u .*= 2)),
        DiscreteCallback((u, t, integrator) -> false, integrator -> nothing)
    )
    wrapper = PropertyForwardingWrapper(nothing, pairs((;)))
    if VERSION >= v"1.12"
        @test get(DiffEqBase.merge_problem_kwargs(wrapper; erase_callbacks = false), :erase_callbacks, nothing) === false
    end
    for run in (init, solve), setting in (false, nothing)
        integ = run(wrapper, Tsit5(); callback = callbacks, erase_callbacks = setting, abstol = 1.0e-10, reltol = 1.0e-10)
        if setting === false
            @test integ.opts.callback.continuous_callbacks isa Tuple
            @test integ.opts.callback.discrete_callbacks isa Tuple
        elseif VERSION >= v"1.12"
            @test integ.opts.callback.continuous_callbacks isa Vector{Any}
        end
        @test solve!(integ).u[end][1] ≈ 2 * exp(-1) atol = 1.0e-8
    end
end
