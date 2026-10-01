using OrdinaryDiffEqSDIRK, Test, SparseArrays, LinearAlgebra

f_oop(u, p, t) = -1000 * u
prob = ODEProblem(f_oop, 1.0, (0.0, 1.0))

sol = solve(prob, TRBDF2(); adaptive = true)
@test sol.retcode == ReturnCode.Success
@test sol.t[end] == 1.0
@test abs(sol.u[end]) < 1.0e-6
@test 0 < sol.stats.njacs < sol.stats.nw

using StaticArrays
f_oop_sa(u, p, t) = SA[-1000.0 * u[1], -500.0 * u[2]]
prob_sa = ODEProblem(f_oop_sa, SA[1.0, 1.0], (0.0, 1.0))

sol_sa = solve(prob_sa, TRBDF2(); adaptive = true)
@test sol_sa.retcode == ReturnCode.Success
@test sol_sa.t[end] == 1.0
@test 0 < sol_sa.stats.njacs < sol_sa.stats.nw

function f_iip!(du, u, p, t)
    du[1] = -1000.0 * u[1]
    du[2] = -500.0 * u[2]
    return nothing
end
prob_iip = ODEProblem(f_iip!, [1.0, 1.0], (0.0, 1.0))
sol_iip = solve(prob_iip, TRBDF2(); adaptive = true)
@test 0 < sol_iip.stats.njacs < sol_iip.stats.nw

A = spdiagm(0 => [-1000.0, -1.0])
M = Diagonal([2.0, 1.0])
f_sparse(u, p, t) = A * u
jac_sparse(u, p, t) = copy(A)
prob_sparse = ODEProblem(
    ODEFunction(f_sparse; jac = jac_sparse, jac_prototype = A, mass_matrix = M),
    [1.0, 1.0], (0.0, 1.0)
)
sol_sparse = solve(prob_sparse, TRBDF2(); reltol = 1.0e-7, abstol = 1.0e-9)
@test sol_sparse.retcode == ReturnCode.Success
@test 0 < sol_sparse.stats.njacs < sol_sparse.stats.nw
@test abs(sol_sparse.u[end][2] - exp(-1)) < 5.0e-5
