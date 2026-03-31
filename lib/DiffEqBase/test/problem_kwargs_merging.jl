using DiffEqBase, Test

# Test merge_problem_kwargs function - only using DiffEqBase directly

@testset "merge_problem_kwargs API" begin
    # Create a simple problem for testing
    function f(du, u, p, t)
        du[1] = -u[1]
    end

    # Test 1: Problem with no kwargs
    prob_no_kwargs = ODEProblem(f, [1.0], (0.0, 1.0))
    kwargs_in = (abstol = 1.0e-6, reltol = 1.0e-6)
    kwargs_out = DiffEqBase.merge_problem_kwargs(prob_no_kwargs; kwargs_in...)
    @test kwargs_out == Base.pairs(kwargs_in)

    # Test 2: Problem with empty kwargs
    prob_empty_kwargs = ODEProblem(f, [1.0], (0.0, 1.0); Dict{Symbol, Any}()...)
    kwargs_out = DiffEqBase.merge_problem_kwargs(prob_empty_kwargs; kwargs_in...)
    @test kwargs_out == Base.pairs(kwargs_in)

    # Test 3: Problem kwargs are preserved, passed kwargs take precedence
    prob_with_kwargs = ODEProblem(f, [1.0], (0.0, 1.0); abstol = 1.0e-8, reltol = 1.0e-8)
    kwargs_out = DiffEqBase.merge_problem_kwargs(prob_with_kwargs; (abstol = 1.0e-6,)...)
    @test kwargs_out.abstol == 1.0e-6  # Passed kwarg takes precedence
    @test kwargs_out.reltol == 1.0e-8  # Problem kwarg is preserved

    # Test 4: Multiple kwargs from problem and passed
    prob_multi = ODEProblem(f, [1.0], (0.0, 1.0); abstol = 1.0e-8, saveat = 0.1)
    kwargs_out = DiffEqBase.merge_problem_kwargs(prob_multi; (reltol = 1.0e-6, maxiters = 1000)...)
    @test kwargs_out.abstol == 1.0e-8    # From problem
    @test kwargs_out.saveat == 0.1     # From problem
    @test kwargs_out.reltol == 1.0e-6    # From passed kwargs
    @test kwargs_out.maxiters == 1000  # From passed kwargs

    # Test 5: Callback merging disabled
    condition(u, t, integrator) = t - 0.5
    affect!(integrator) = (integrator.u[1] += 1.0)
    cb1 = ContinuousCallback(condition, affect!)
    cb2 = ContinuousCallback(condition, affect!)

    prob_with_cb = ODEProblem(f, [1.0], (0.0, 1.0); callback = cb1)
    kwargs_out = DiffEqBase.merge_problem_kwargs(
        prob_with_cb; merge_callbacks = false,
        (callback = cb2,)...
    )
    @test kwargs_out.callback === cb2  # cb2 should override cb1

    # Test 6: Callback merging enabled (default)
    kwargs_out = DiffEqBase.merge_problem_kwargs(prob_with_cb; (callback = cb2,)...)
    @test kwargs_out.callback isa DiffEqBase.CallbackSet
    @test length(kwargs_out.callback.continuous_callbacks) == 2

    # Test 7: Callback merging with CallbackSet
    cb3 = ContinuousCallback(condition, affect!)
    cbset = CallbackSet(cb3)
    kwargs_out = DiffEqBase.merge_problem_kwargs(prob_with_cb; (callback = cbset,)...)
    @test kwargs_out.callback isa DiffEqBase.CallbackSet
    @test length(kwargs_out.callback.continuous_callbacks) == 2

    # Test 8: Only problem has callback
    kwargs_out = DiffEqBase.merge_problem_kwargs(prob_with_cb; (abstol = 1.0e-6,)...)
    @test kwargs_out.callback === cb1
    @test kwargs_out.abstol == 1.0e-6

    # Test 9: Only passed kwargs have callback
    prob_no_cb = ODEProblem(f, [1.0], (0.0, 1.0); abstol = 1.0e-8)
    kwargs_out = DiffEqBase.merge_problem_kwargs(prob_no_cb; (callback = cb2, reltol = 1.0e-6)...)
    @test kwargs_out.callback === cb2
    @test kwargs_out.abstol == 1.0e-8
    @test kwargs_out.reltol == 1.0e-6

    # Test 10: DiscreteCallback merging
    dcb1 = DiscreteCallback((u, t, integrator) -> true, integrator -> nothing)
    dcb2 = DiscreteCallback((u, t, integrator) -> true, integrator -> nothing)
    prob_with_dcb = ODEProblem(f, [1.0], (0.0, 1.0); callback = dcb1)
    kwargs_out = DiffEqBase.merge_problem_kwargs(prob_with_dcb; (callback = dcb2,)...)
    @test kwargs_out.callback isa DiffEqBase.CallbackSet
    @test length(kwargs_out.callback.discrete_callbacks) == 2

    # Test 11: Mixed Continuous and Discrete callback merging
    prob_with_ccb = ODEProblem(f, [1.0], (0.0, 1.0); callback = cb1)
    kwargs_out = DiffEqBase.merge_problem_kwargs(prob_with_ccb; (callback = dcb1,)...)
    @test kwargs_out.callback isa DiffEqBase.CallbackSet
    @test length(kwargs_out.callback.continuous_callbacks) == 1
    @test length(kwargs_out.callback.discrete_callbacks) == 1
end
