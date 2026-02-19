using OrdinaryDiffEqBDF, DiffEqBase, Test

# Test that u_modified! triggers DAE initialization
# Regression test for https://github.com/SciML/DifferentialEquations.jl/issues/1127

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
