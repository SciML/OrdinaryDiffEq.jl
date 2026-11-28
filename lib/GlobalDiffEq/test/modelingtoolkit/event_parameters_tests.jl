using GlobalDiffEq, ModelingToolkit, OrdinaryDiffEqTsit5, SciMLSensitivity, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

relerr(estimate, actual) = abs(estimate - actual) / abs(actual)

# x' = -k x with a discrete parameter k that an event sets from the state. The parameters
# are an `MTKParameters` (in a `DespecializedParameters`), so GlobalAdjoint reads and sets
# the event parameter with SymbolicIndexingInterface (`getp`, `remake_buffer`).
@variables x(t) = 1.0
@discretes k(t) = 1.0
τ = log(2.0)
cases = (
    # At t = 1, k = -3x. Exact x(2) = exp(-1) exp(3 exp(-1)).
    (
        ModelingToolkit.SymbolicDiscreteCallback([1.0] => [k ~ -3 * Pre(x)]; discrete_parameters = [k]),
        :discrete_events, exp(-1) * exp(3 * exp(-1)),
    ),
    # At x = 0.5 (t = log 2), k = 4x = 2. Exact x(2) = 0.5 exp(-2 (2 - log 2)).
    (
        ModelingToolkit.SymbolicContinuousCallback([x ~ 0.5] => [k ~ 4 * Pre(x)]; discrete_parameters = [k]),
        :continuous_events, 0.5 * exp(-2 * (2 - τ)),
    ),
)

@testset "GlobalAdjoint: an MTK discrete parameter set from the state ($kind)" for (event, kind, exact) in cases
    events = kind === :discrete_events ? (; discrete_events = [event]) : (; continuous_events = [event])
    @named model = System([D(x) ~ -k * x], t, [x, k], []; events...)
    prob = ODEProblem(mtkcompile(model), [], (0.0, 2.0))
    for tol in (1.0e-4, 1.0e-6)
        sol = solve(deepcopy(prob), Tsit5(); abstol = tol, reltol = tol)
        actual = abs(sol[x][end] - exact)
        # The VJP products go through `remake_buffer` without a fallback warning.
        vjp = @test_logs min_level = Base.CoreLogging.Warn adjoint_error_estimate(
            deepcopy(prob), GlobalAdjoint(Tsit5()); abstol = tol, reltol = tol
        )
        finite = adjoint_error_estimate(
            deepcopy(prob), GlobalAdjoint(Tsit5(); sensealg = GaussKronrodAdjoint(autojacvec = false));
            abstol = tol, reltol = tol
        )
        @test relerr(vjp, actual) < 1.0e-3
        @test relerr(finite, actual) < 1.0e-3
    end
end
