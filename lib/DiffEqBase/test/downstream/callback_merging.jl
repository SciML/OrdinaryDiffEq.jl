using OrdinaryDiffEq, DiffEqBase, SciMLBase, Test
using OrdinaryDiffEqLowOrderRK: Euler

# Basic auto callback merging test
do_nothing = DiscreteCallback(
    (u, t, integrator) -> true,
    integrator -> nothing
)
problem = ODEProblem(
    (u, p, t) -> -u,
    1.0, (0.0, 1.0),
    callback = do_nothing
)
solve(
    problem, Euler(),
    dt = 0.1,
    callback = do_nothing
)

@testset "Callback merging through solve" begin
    # Test that callbacks passed to problem constructor are properly merged
    # with callbacks passed to solve

    function lorenz!(du, u, p, t)
        du[1] = 10.0(u[2] - u[1])
        du[2] = u[1] * (28.0 - u[3]) - u[2]
        du[3] = u[1] * u[2] - (8 / 3) * u[3]
    end
    u0 = [1.0, 0.0, 0.0]
    tspan = (0.0, 1.0)

    # Create two callbacks that track when they're called
    cb1_called = Ref(false)
    condition1(u, t, integrator) = t - 0.3
    affect1!(integrator) = (cb1_called[] = true)
    cb1 = ContinuousCallback(condition1, affect1!)

    cb2_called = Ref(false)
    condition2(u, t, integrator) = t - 0.7
    affect2!(integrator) = (cb2_called[] = true)
    cb2 = ContinuousCallback(condition2, affect2!)

    # Test 1: Callback in problem constructor only
    cb1_called[] = false
    prob = ODEProblem(lorenz!, u0, tspan; callback = cb1)
    sol = solve(prob, Tsit5())
    @test cb1_called[]
    @test sol.t[end] ≈ 1.0

    # Test 2: Callback in solve only
    cb2_called[] = false
    prob = ODEProblem(lorenz!, u0, tspan)
    sol = solve(prob, Tsit5(); callback = cb2)
    @test cb2_called[]
    @test sol.t[end] ≈ 1.0

    # Test 3: Callbacks in both (should merge by default)
    cb1_called[] = false
    cb2_called[] = false
    prob = ODEProblem(lorenz!, u0, tspan; callback = cb1)
    sol = solve(prob, Tsit5(); callback = cb2)
    @test cb1_called[]
    @test cb2_called[]
    @test sol.t[end] ≈ 1.0

    # Test 4: merge_callbacks = false (solve callback should override)
    cb1_called[] = false
    cb2_called[] = false
    prob = ODEProblem(lorenz!, u0, tspan; callback = cb1)
    sol = solve(prob, Tsit5(); callback = cb2, merge_callbacks = false)
    @test !cb1_called[]  # cb1 should not be called
    @test cb2_called[]   # cb2 should be called
    @test sol.t[end] ≈ 1.0
end

@testset "Callback merging through init" begin
    # Test that callbacks are properly merged when using init instead of solve

    function simple!(du, u, p, t)
        du[1] = -u[1]
    end
    u0 = [1.0]
    tspan = (0.0, 1.0)

    cb1_called = Ref(false)
    condition1(u, t, integrator) = t - 0.3
    affect1!(integrator) = (cb1_called[] = true)
    cb1 = ContinuousCallback(condition1, affect1!)

    cb2_called = Ref(false)
    condition2(u, t, integrator) = t - 0.7
    affect2!(integrator) = (cb2_called[] = true)
    cb2 = ContinuousCallback(condition2, affect2!)

    # Test: Callbacks in both problem and init (should merge)
    cb1_called[] = false
    cb2_called[] = false
    prob = ODEProblem(simple!, u0, tspan; callback = cb1)
    integrator = init(prob, Tsit5(); callback = cb2)
    solve!(integrator)
    @test cb1_called[]
    @test cb2_called[]
end

@testset "Other kwargs merging" begin
    # Test that non-callback kwargs are properly merged by checking they're accessible

    function simple!(du, u, p, t)
        du[1] = -u[1]
    end
    u0 = [1.0]
    tspan = (0.0, 1.0)

    # Test that problem kwargs are preserved and solve kwargs override
    prob = ODEProblem(simple!, u0, tspan; abstol = 1.0e-10, saveat = 0.1)

    # Both abstol and saveat should be used
    sol = solve(prob, Tsit5())
    @test sol.t ≈ [i * 0.1 for i in 0:10]  # saveat from problem kwargs

    # Override saveat, keep abstol
    sol = solve(prob, Tsit5(); saveat = 0.5)
    @test sol.t ≈ [0.0, 0.5, 1.0]  # saveat from solve kwargs

    # Test with save_everystep
    prob_everystep = ODEProblem(simple!, u0, tspan; save_everystep = false, saveat = 0.5)
    sol = solve(prob_everystep, Tsit5())
    @test sol.t ≈ [0.0, 0.5, 1.0]  # Only saved at saveat times, not every step

    # Override save_everystep
    sol = solve(prob_everystep, Tsit5(); save_everystep = true)
    @test length(sol.t) > 3  # Should save at more than just saveat times
end

@testset "kwargshandle merging" begin
    # Test that kwargshandle is properly respected from problem kwargs

    function simple!(du, u, p, t)
        du[1] = -u[1]
    end
    u0 = [1.0]
    tspan = (0.0, 1.0)

    # Problem with KeywordArgSilent should not warn on invalid kwargs
    prob = ODEProblem(
        simple!, u0, tspan;
        kwargshandle = SciMLBase.KeywordArgSilent,
        invalid_kwarg = "should be ignored"
    )
    @test_nowarn sol = solve(prob, Tsit5())

    # Problem with KeywordArgWarn and invalid kwarg should warn
    prob = ODEProblem(
        simple!, u0, tspan;
        kwargshandle = SciMLBase.KeywordArgWarn,
        invalid_kwarg = "should warn"
    )
    @test_logs (:warn, SciMLBase.KWARGWARN_MESSAGE) sol = solve(prob, Tsit5())

    # Default should error on invalid kwargs
    prob = ODEProblem(simple!, u0, tspan; invalid_kwarg = "should error")
    @test_throws SciMLBase.CommonKwargError sol = solve(prob, Tsit5())
end

@testset "derivative_discontinuity! sticky-on-true across callbacks" begin
    # When several callbacks fire on the same step, a derivative discontinuity
    # flagged by *any* callback must survive, regardless of evaluation order and
    # regardless of other callbacks declaring no discontinuity. The integrator
    # field is sticky-on-true within a step; an explicit `false` only sticks when
    # no callback flagged `true`. If no callback calls the setter at all, the
    # safe default (treat as discontinuity) is preserved.
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))

    flag_true = DiscreteCallback(
        (u, t, integ) -> true,
        integ -> derivative_discontinuity!(integ, true)
    )
    flag_false = DiscreteCallback(
        (u, t, integ) -> true,
        integ -> derivative_discontinuity!(integ, false)
    )
    silent = DiscreteCallback((u, t, integ) -> true, integ -> nothing)

    # Helper: step once and read the post-callback step-level state.
    function step_flags(cbset)
        integ = init(prob, Tsit5(); callback = cbset)
        step!(integ)
        return integ.derivative_discontinuity, integ.user_set_discontinuity
    end

    # true wins regardless of order
    @test step_flags(CallbackSet(flag_true, flag_false)) == (true, true)
    @test step_flags(CallbackSet(flag_false, flag_true)) == (true, true)

    # all-false: explicit false sticks, no recomputation flagged
    @test step_flags(CallbackSet(flag_false, flag_false)) == (false, true)

    # no setter call at all: safe default keeps discontinuity, marked unset
    @test step_flags(CallbackSet(silent)) == (true, false)

    # a silent callback does not suppress another callback's explicit false
    @test step_flags(CallbackSet(silent, flag_false)) == (true, true)

    # full solves complete in both orders
    @test solve(prob, Tsit5(); callback = CallbackSet(flag_true, flag_false)).retcode ==
        ReturnCode.Success
    @test solve(prob, Tsit5(); callback = CallbackSet(flag_false, flag_true)).retcode ==
        ReturnCode.Success
end
