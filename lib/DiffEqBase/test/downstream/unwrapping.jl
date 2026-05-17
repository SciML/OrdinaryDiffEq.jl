using OrdinaryDiffEq, OrdinaryDiffEqSDIRK, SciMLBase, Test
using FunctionWrappersWrappers
my_f(u, p, t) = u
my_f!(du, u, p, t) = du .= u
ode = ODEProblem(my_f, [1.0], (0.0, 1.0))
integrator = init(ode, Tsit5())
@test SciMLBase.unwrapped_f(integrator.f.f) === my_f

ode = ODEProblem(my_f!, [1.0], (0.0, 1.0))
integrator = init(ode, Tsit5())
@test SciMLBase.unwrapped_f(integrator.f.f) === my_f!

using OrdinaryDiffEq, ForwardDiff, Measurements
x = 1.0 ± 0.0
f = (du, u, p, t) -> du .= u
tspan = (0.0, 1.0)
prob = ODEProblem(f, [x], tspan)

# Should not error during problem construction but should be unwrapped
integ = init(prob, Tsit5(), dt = 0.1)
@test SciMLBase.unwrapped_f(integ.f.f) === f

# Handle functional initial conditions
prob = ODEProblem((dx, x, p, t) -> (dx .= 0), (p, t) -> zeros(2), (0, 10))
solve(prob, TRBDF2())

# AutoSpecialize skips the FunctionWrappersWrapper when u0's eltype is already
# a ForwardDiff.Dual. The wrapper's signature slots are compiled for the seen
# input types; under an outer ForwardDiff layer the inner solver re-widens
# `u`/`p` into deeper-nested Duals that the wrapped slot doesn't match,
# producing spurious `NoFunctionWrapperFoundError` or dispatch into a slot
# whose compiled body produces wrong-typed Duals. Skipping the wrap lets
# Julia specialize directly on the actual nested-Dual types.
let
    g!(du, u, p, t) = (du .= -u; nothing)
    Tag = ForwardDiff.Tag{:UnwrappingDualU0Test, Float64}
    d(x) = ForwardDiff.Dual{Tag, Float64, 1}(x, ForwardDiff.Partials{1, Float64}((1.0,)))

    ode_auto = ODEFunction{true, SciMLBase.AutoSpecialize}(g!)
    prob_dual = ODEProblem(ode_auto, [d(1.0), d(1.0)], (0.0, 1.0), [d(1.0), d(1.0)])
    integ_dual = init(prob_dual, Tsit5())
    @test SciMLBase.unwrapped_f(integ_dual.f.f) === g!
    @test integ_dual.f.f === g!

    # Regression guard: Float64 u0 still goes through the wrapper.
    prob_float = ODEProblem(ode_auto, [1.0, 1.0], (0.0, 1.0), [1.0, 1.0])
    integ_float = init(prob_float, Tsit5())
    @test integ_float.f.f isa
        FunctionWrappersWrappers.FunctionWrappersWrapper
end
