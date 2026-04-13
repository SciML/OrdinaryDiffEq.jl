using OrdinaryDiffEqBDF, DiffEqBase, Test
import OrdinaryDiffEqCore: _initialize_dae!

# Test that u_modified! triggers DAE initialization
# Regression test for https://github.com/SciML/DifferentialEquations.jl/issues/1127

# Custom init algorithm that counts how many times initialize_dae! is called
struct CountingInit{T} <: SciMLBase.DAEInitializationAlgorithm
    alg::T
    count::Ref{Int}
end

function _initialize_dae!(integrator, prob, alg::CountingInit, isinplace)
    alg.count[] += 1
    return _initialize_dae!(integrator, prob, alg.alg, isinplace)
end

function f(out, du, u, _, _)
    out[1] = -0.04u[1] + 1.0e4 * u[2] * u[3] - du[1]
    out[2] = +0.04u[1] - 3.0e7 * u[2]^2 - 1.0e4 * u[2] * u[3] - du[2]
    return out[3] = u[1] + u[2] + u[3] - 1.0
end

u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]

prob = DAEProblem(f, du₀, u₀, (0.0, 1.0), differential_vars = [true, true, false])

# With CheckInit, modifying u to break constraints and calling u_modified!
# should trigger the initialization check and error
int = init(prob, DFBDF(), initializealg = DiffEqBase.CheckInit())
int.u[1] = 2.0 # Breaks algebraic constraint u[1] + u[2] + u[3] = 1
u_modified!(int, true)
@test_throws SciMLBase.CheckInitFailureError step!(int)

# With default init, modifying u and calling u_modified! should
# reinitialize the algebraic variables to satisfy constraints
int2 = init(prob, DFBDF())
int2.u[1] = 2.0 # Breaks algebraic constraint
u_modified!(int2, true)
step!(int2)
@test int2.u[1] + int2.u[2] + int2.u[3] ≈ 1.0 atol = 1.0e-10

# Test that u_modified! calls initialize_dae! exactly once
counter = CountingInit(DiffEqBase.CheckInit(), Ref(0))
int3 = init(prob, DFBDF(), initializealg = counter)
init_after_init = counter.count[]  # init() calls initialize_dae! once
@test init_after_init == 1
# Modify u without breaking constraints (just set same values)
int3.u .= u₀
u_modified!(int3, true)
step!(int3)
@test counter.count[] == init_after_init + 1  # exactly one more init call

# Test that callbacks don't trigger double initialization.
# A callback that modifies u should trigger initialize_dae! exactly once
# (via reeval_internals_due_to_modification!, NOT again in loopheader!)
prob2 = DAEProblem(
    f, du₀, u₀, (0.0, 100.0),
    differential_vars = [true, true, false]
)
counter2 = CountingInit(DiffEqBase.CheckInit(), Ref(0))
condition(u, t, integrator) = true  # fire every step
affect!(integrator) = (integrator.u .= integrator.u)  # no-op modification
cb = DiscreteCallback(condition, affect!)
int4 = init(prob2, DFBDF(), initializealg = counter2, callback = cb)
init_count_before = counter2.count[]
step!(int4)
# Callback fires once → reeval_internals_due_to_modification! → initialize_dae! once
# If double-init existed, this would be 2
@test counter2.count[] == init_count_before + 1
