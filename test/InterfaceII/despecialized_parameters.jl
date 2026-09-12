using OrdinaryDiffEq, Test
import DiffEqBase
import SciMLBase

@test OrdinaryDiffEq.AutoDespecialize === SciMLBase.AutoDespecialize
@test OrdinaryDiffEq.AutoRespecialize === SciMLBase.AutoRespecialize
@test OrdinaryDiffEq.AutoDePSpecialize === SciMLBase.AutoDePSpecialize

struct FirstParameterLayout
    rate::Float64
end

struct SecondParameterLayout
    rate::Float64
    unused::Int
end

function despecialized_decay!(du, u, p, t)
    du[1] = -p.rate * u[1]
    return
end

function despecialized_problem(p)
    return ODEProblem{true, SciMLBase.AutoDespecialize}(
        despecialized_decay!, [1.0], (0.0, 1.0), p
    )
end

first_prob = despecialized_problem(FirstParameterLayout(1.0))
second_prob = despecialized_problem(SecondParameterLayout(2.0, 1))

@test first_prob.p isa FirstParameterLayout
@test second_prob.p isa SecondParameterLayout

first_concrete = DiffEqBase.get_concrete_problem(first_prob, true; alg = Tsit5())
second_concrete = DiffEqBase.get_concrete_problem(second_prob, true; alg = Tsit5())
@test typeof(first_concrete.p) === typeof(second_concrete.p) ===
    SciMLBase.DespecializedParameters
@test typeof(first_concrete) === typeof(second_concrete)
@test solve(
    first_prob, Tsit5(); saveat = [1.0], abstol = 1.0e-10, reltol = 1.0e-10
).u[end][1] ≈ exp(-1)
@test solve(
    second_prob, Tsit5(); saveat = [1.0], abstol = 1.0e-10, reltol = 1.0e-10
).u[end][1] ≈ exp(-2)
@test solve(
    first_prob, Rosenbrock23(); saveat = [1.0], abstol = 1.0e-12, reltol = 1.0e-12
).u[end][1] ≈ exp(-1)

function callback_constant_rhs!(du, u, p, t)
    fill!(du, 0)
    return nothing
end

function method_instance_count(method)
    return count(instance -> instance isa Core.MethodInstance, Base.specializations(method))
end

Base.@noinline function solve_without_specializing!(integrator)
    Base.@nospecialize integrator
    # Keep the test harness from compiling `solve!` before its method instances are counted.
    return Base.invokelatest(SciMLBase.solve!, integrator)
end

@testset "callback changes do not compile another solver" begin
    first_callback = DiscreteCallback(
        (u, t, integrator) -> iszero(u[1]),
        integrator -> (integrator.u[1] = 1.0)
    )
    second_callback = DiscreteCallback(
        (u, t, integrator) -> iszero(u[1]),
        integrator -> (integrator.u[1] = 2.0)
    )

    for specialize in (
            SciMLBase.AutoSpecialize,
            SciMLBase.AutoDespecialize,
            SciMLBase.NoSpecialize,
        )
        first_problem = ODEProblem{true, specialize}(
            callback_constant_rhs!, [0.0], (0.0, 1.0); callback = first_callback
        )
        second_problem = ODEProblem{true, specialize}(
            callback_constant_rhs!, [0.0], (0.0, 1.0); callback = second_callback
        )
        first_integrator = init(first_problem, Tsit5())
        second_integrator = init(second_problem, Tsit5())
        solve_method = which(SciMLBase.solve!, Tuple{typeof(first_integrator)})

        solve_without_specializing!(first_integrator)
        first_count = method_instance_count(solve_method)
        solve_without_specializing!(second_integrator)
        @test method_instance_count(solve_method) == first_count
    end
end

@testset "callback types do not specialize despecialized integrators" begin
    first_callback = DiscreteCallback(
        (u, t, integrator) -> iszero(u[1]),
        integrator -> (integrator.u[1] = 1.0)
    )
    second_callback = DiscreteCallback(
        (u, t, integrator) -> iszero(u[1]),
        integrator -> (integrator.u[1] = 2.0)
    )
    unused_callback = DiscreteCallback(
        (u, t, integrator) -> false,
        integrator -> nothing
    )

    for specialize in (SciMLBase.AutoDespecialize, SciMLBase.NoSpecialize)
        no_callback_problem = ODEProblem{true, specialize}(
            callback_constant_rhs!, [0.0], (0.0, 1.0)
        )
        first_callback_problem = ODEProblem{true, specialize}(
            callback_constant_rhs!, [0.0], (0.0, 1.0); callback = first_callback
        )
        second_callback_problem = ODEProblem{true, specialize}(
            callback_constant_rhs!, [0.0], (0.0, 1.0);
            callback = CallbackSet(second_callback, unused_callback)
        )

        no_callback_concrete = DiffEqBase.get_concrete_problem(
            no_callback_problem, true; alg = Tsit5()
        )
        first_concrete = DiffEqBase.get_concrete_problem(
            first_callback_problem, true; alg = Tsit5()
        )
        second_concrete = DiffEqBase.get_concrete_problem(
            second_callback_problem, true; alg = Tsit5()
        )
        @test first_concrete.kwargs[:callback].continuous_callbacks isa Vector{Any}
        @test first_concrete.kwargs[:callback].discrete_callbacks isa Vector{Any}
        @test typeof(no_callback_concrete) === typeof(first_concrete)
        @test typeof(first_concrete) === typeof(second_concrete)

        no_callback_integrator = init(no_callback_problem, Tsit5())
        first_integrator = init(first_callback_problem, Tsit5())
        second_integrator = init(second_callback_problem, Tsit5())
        first_solve_callback_integrator = init(
            no_callback_problem, Tsit5(); callback = first_callback
        )
        second_solve_callback_integrator = init(
            no_callback_problem, Tsit5(); callback = CallbackSet(second_callback, unused_callback)
        )
        @test typeof(no_callback_integrator) === typeof(first_integrator)
        @test typeof(first_integrator) === typeof(second_integrator)
        @test first_solve_callback_integrator.opts.callback.discrete_callbacks isa Vector{Any}
        @test typeof(first_integrator) === typeof(first_solve_callback_integrator)
        @test typeof(first_solve_callback_integrator) === typeof(second_solve_callback_integrator)
        @test solve(first_callback_problem, Tsit5()).u[end] == [1.0]
        @test solve(second_callback_problem, Tsit5()).u[end] == [2.0]
        @test solve(no_callback_problem, Tsit5(); callback = first_callback).u[end] == [1.0]

        first_continuous = ContinuousCallback(
            (u, t, integrator) -> t - 0.5,
            integrator -> (integrator.u[1] += 1)
        )
        early_continuous = ContinuousCallback(
            (u, t, integrator) -> t - 0.25,
            integrator -> (integrator.u[1] += 1)
        )
        late_continuous = ContinuousCallback(
            (u, t, integrator) -> t - 0.75,
            integrator -> (integrator.u[1] += 1)
        )
        first_continuous_problem = ODEProblem{true, specialize}(
            callback_constant_rhs!, [0.0], (0.0, 1.0); callback = first_continuous
        )
        second_continuous_problem = ODEProblem{true, specialize}(
            callback_constant_rhs!, [0.0], (0.0, 1.0);
            callback = CallbackSet(early_continuous, late_continuous)
        )

        first_continuous_integrator = init(first_continuous_problem, Tsit5())
        second_continuous_integrator = init(second_continuous_problem, Tsit5())
        @test typeof(first_continuous_integrator) === typeof(second_continuous_integrator)
        @test solve(first_continuous_problem, Tsit5()).u[end] == [1.0]
        @test solve(second_continuous_problem, Tsit5()).u[end] == [2.0]

        first_vector_continuous = VectorContinuousCallback(
            (out, u, t, integrator) -> (out[1] = t - 0.5),
            (integrator, event_index) -> (integrator.u[1] += 1),
            1
        )
        early_vector_continuous = VectorContinuousCallback(
            (out, u, t, integrator) -> (out[1] = t - 0.25),
            (integrator, event_index) -> (integrator.u[1] += 1),
            1
        )
        late_vector_continuous = VectorContinuousCallback(
            (out, u, t, integrator) -> (out[1] = t - 0.75),
            (integrator, event_index) -> (integrator.u[1] += 1),
            1
        )
        first_vector_problem = ODEProblem{true, specialize}(
            callback_constant_rhs!, [0.0], (0.0, 1.0); callback = first_vector_continuous
        )
        second_vector_problem = ODEProblem{true, specialize}(
            callback_constant_rhs!, [0.0], (0.0, 1.0);
            callback = CallbackSet(early_vector_continuous, late_vector_continuous)
        )

        first_vector_integrator = init(first_vector_problem, Tsit5())
        second_vector_integrator = init(second_vector_problem, Tsit5())
        @test typeof(first_vector_integrator) === typeof(second_vector_integrator)
        @test solve(first_vector_problem, Tsit5()).u[end] == [1.0]
        @test solve(second_vector_problem, Tsit5()).u[end] == [2.0]
    end
end

@testset "AutoSpecialize wraps callbacks after constructing the integrator" begin
    first_discrete = DiscreteCallback(
        (u, t, integrator) -> iszero(u[1]),
        integrator -> (integrator.u[1] = 1.0)
    )
    second_discrete = DiscreteCallback(
        (u, t, integrator) -> iszero(u[1]),
        integrator -> (integrator.u[1] = 2.0)
    )
    unused_discrete = DiscreteCallback(
        (u, t, integrator) -> false,
        integrator -> nothing
    )
    no_callback_problem = ODEProblem{true, SciMLBase.AutoSpecialize}(
        callback_constant_rhs!, [0.0], (0.0, 1.0)
    )
    first_discrete_problem = ODEProblem{true, SciMLBase.AutoSpecialize}(
        callback_constant_rhs!, [0.0], (0.0, 1.0); callback = first_discrete
    )
    second_discrete_problem = ODEProblem{true, SciMLBase.AutoSpecialize}(
        callback_constant_rhs!, [0.0], (0.0, 1.0);
        callback = CallbackSet(second_discrete, unused_discrete)
    )

    no_callback_integrator = init(no_callback_problem, Tsit5())
    first_discrete_integrator = init(first_discrete_problem, Tsit5())
    second_discrete_integrator = init(second_discrete_problem, Tsit5())
    first_wrapped = only(first_discrete_integrator.opts.callback.discrete_callbacks)
    second_wrapped = first(second_discrete_integrator.opts.callback.discrete_callbacks)
    @test typeof(no_callback_integrator) === typeof(first_discrete_integrator)
    @test typeof(first_discrete_integrator) === typeof(second_discrete_integrator)
    @test typeof(first_wrapped.condition) === typeof(second_wrapped.condition)
    @test typeof(first_wrapped.affect!) === typeof(second_wrapped.affect!)
    @test solve(first_discrete_problem, Tsit5()).u[end] == [1.0]
    @test solve(second_discrete_problem, Tsit5()).u[end] == [2.0]

    first_continuous = ContinuousCallback(
        (u, t, integrator) -> t - 0.5,
        integrator -> (integrator.u[1] += 1)
    )
    second_continuous = ContinuousCallback(
        (u, t, integrator) -> t - 0.25,
        integrator -> (integrator.u[1] += 1)
    )
    first_continuous_integrator = init(
        no_callback_problem, Tsit5(); callback = first_continuous
    )
    second_continuous_integrator = init(
        no_callback_problem, Tsit5(); callback = second_continuous
    )
    first_continuous_wrapped = only(
        first_continuous_integrator.opts.callback.continuous_callbacks
    )
    second_continuous_wrapped = only(
        second_continuous_integrator.opts.callback.continuous_callbacks
    )
    @test typeof(first_continuous_integrator) === typeof(second_continuous_integrator)
    @test typeof(first_continuous_wrapped.condition) ===
        typeof(second_continuous_wrapped.condition)
    @test typeof(first_continuous_wrapped.affect!) ===
        typeof(second_continuous_wrapped.affect!)
    @test solve(no_callback_problem, Tsit5(); callback = first_continuous).u[end] == [1.0]

    negative_only = ContinuousCallback(
        (u, t, integrator) -> 0.5 - t,
        nothing,
        integrator -> (integrator.u[1] += 1)
    )
    negative_only_integrator = init(no_callback_problem, Tsit5(); callback = negative_only)
    negative_only_wrapped = only(
        negative_only_integrator.opts.callback.continuous_callbacks
    )
    @test negative_only_wrapped.affect! === nothing
    @test solve(no_callback_problem, Tsit5(); callback = negative_only).u[end] == [1.0]

    first_vector = VectorContinuousCallback(
        (out, u, t, integrator) -> (out[1] = t - 0.5),
        (integrator, events) -> (integrator.u[1] += count(!iszero, events)),
        1
    )
    second_vector = VectorContinuousCallback(
        (out, u, t, integrator) -> (out[1] = t - 0.25),
        (integrator, events) -> (integrator.u[1] += count(!iszero, events)),
        1
    )
    first_vector_integrator = init(no_callback_problem, Tsit5(); callback = first_vector)
    second_vector_integrator = init(no_callback_problem, Tsit5(); callback = second_vector)
    first_vector_wrapped = only(first_vector_integrator.opts.callback.continuous_callbacks)
    second_vector_wrapped = only(second_vector_integrator.opts.callback.continuous_callbacks)
    @test typeof(first_vector_integrator) === typeof(second_vector_integrator)
    @test typeof(first_vector_wrapped.condition) === typeof(second_vector_wrapped.condition)
    @test typeof(first_vector_wrapped.affect!) === typeof(second_vector_wrapped.affect!)
    @test solve(no_callback_problem, Tsit5(); callback = first_vector).u[end] == [1.0]

    no_vector_affect = VectorContinuousCallback(
        (out, u, t, integrator) -> (out[1] = t - 0.5), nothing, 1
    )
    no_vector_affect_integrator = init(
        no_callback_problem, Tsit5(); callback = no_vector_affect
    )
    no_vector_affect_wrapped = only(
        no_vector_affect_integrator.opts.callback.continuous_callbacks
    )
    @test no_vector_affect_wrapped.affect! === nothing
    @test solve(no_callback_problem, Tsit5(); callback = no_vector_affect).u[end] == [0.0]
end
