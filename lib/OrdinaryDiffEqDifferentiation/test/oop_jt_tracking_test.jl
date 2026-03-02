# Test that OOP update_W! correctly tracks J_t when new_W=true but new_jac=false.
# Without the fix, OOP TRBDF2 can hang in an infinite @goto REDO loop because
# isJcurrent() returns false even after J was freshly computed by calc_W.
using OrdinaryDiffEqSDIRK, Test

# Stiff OOP scalar problem: du/dt = -1000*u
# The extreme stiffness forces step retries where new_W=true, new_jac=false can occur.
f_oop(u, p, t) = -1000 * u
prob = ODEProblem(f_oop, 1.0, (0.0, 1.0))

# This should complete without hanging. The bug causes an infinite loop in nlsolve.
sol = solve(prob, TRBDF2(); adaptive = true)
@test sol.retcode == ReturnCode.Success
@test sol.t[end] == 1.0
@test abs(sol.u[end]) < 1.0e-6  # Should decay to near zero

# Also test with a vector OOP problem (StaticArrays-style)
using StaticArrays
f_oop_sa(u, p, t) = SA[-1000.0 * u[1], -500.0 * u[2]]
prob_sa = ODEProblem(f_oop_sa, SA[1.0, 1.0], (0.0, 1.0))

sol_sa = solve(prob_sa, TRBDF2(); adaptive = true)
@test sol_sa.retcode == ReturnCode.Success
@test sol_sa.t[end] == 1.0
