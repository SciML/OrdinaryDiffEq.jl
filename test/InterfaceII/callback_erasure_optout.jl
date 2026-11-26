using OrdinaryDiffEqTsit5, SciMLBase, Test

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
