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

vdp_oop(u, μ, t) = [μ * ((1 - u[2]^2) * u[1] - u[2]), u[1]]
vdp_iip!(du, u, μ, t) = (du .= vdp_oop(u, μ, t); nothing)
vdp_u0 = [sqrt(3.0), 0.0]
vdp_oop_sol = solve(
    ODEProblem(vdp_oop, vdp_u0, (0.0, 1.0), 1.0e6), ImplicitEuler();
    reltol = 1.0e-6, abstol = 1.0e-8
)
vdp_iip_sol = solve(
    ODEProblem(vdp_iip!, copy(vdp_u0), (0.0, 1.0), 1.0e6), ImplicitEuler();
    reltol = 1.0e-6, abstol = 1.0e-8
)
@test vdp_oop_sol.retcode == vdp_iip_sol.retcode == ReturnCode.Success
@test vdp_oop_sol.stats.njacs == vdp_iip_sol.stats.njacs
@test vdp_oop_sol.stats.nw == vdp_iip_sol.stats.nw
@test vdp_oop_sol.stats.naccept == vdp_iip_sol.stats.naccept
@test vdp_oop_sol.stats.nreject == vdp_iip_sol.stats.nreject

function bertolazzi_rhs(u, p, t)
    f1 = 5 * u[2] * u[3] / (1.0e-2 + (u[2] * u[3])^2) +
        u[2] * u[3] / (1.0e-16 + u[2] * u[3] * (1.0e-8 + u[2] * u[3]))
    f2 = 10 * u[1] * u[3]^2
    f3 = 0.1 * (u[3] - u[2] - 2.5)^2 * u[1] * u[2]
    return SA[2 * f1 - f2 - f3, -f1 + 2 * f2 - f3, -f1 - f2 + 2 * f3]
end

bertolazzi_vec(u, p, t) = Vector(bertolazzi_rhs(SVector{3}(u), p, t))
bertolazzi_iip!(du, u, p, t) = (du .= bertolazzi_rhs(SVector{3}(u), p, t); nothing)

@testset "OOP conservation with Jacobian reuse" begin
    u0 = [0.0, 1.0, 2.0]
    svprob = ODEProblem(bertolazzi_rhs, SVector{3}(u0), (0.0, 1.0))
    vecprob = ODEProblem(bertolazzi_vec, u0, (0.0, 1.0))
    iipprob = ODEProblem(bertolazzi_iip!, copy(u0), (0.0, 1.0))
    for prob in (svprob, vecprob, iipprob)
        sol = solve(prob, TRBDF2(); dt = 1.0e-3)
        @test sol.retcode == ReturnCode.Success
        @test abs(sum(sol.u[end]) - 3.0) < 1.0e-6
    end
    sol = solve(svprob, ImplicitEuler(); dt = 1.0e-3)
    @test sol.retcode == ReturnCode.Success
    @test abs(sum(sol.u[end]) - 3.0) < 1.0e-6
end

function long_rober!(du, u, p, t)
    du[1] = -0.04u[1] + 1.0e4u[2] * u[3]
    du[2] = 0.04u[1] - 1.0e4u[2] * u[3] - 3.0e7u[2]^2
    du[3] = 3.0e7u[2]^2
    return nothing
end

@testset "in-place long-horizon Jacobian reuse" begin
    prob = ODEProblem(long_rober!, [1.0, 0.0, 0.0], (0.0, 1.0e11))
    for alg in (TRBDF2(), KenCarp4())
        sol = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-10)
        @test sol.retcode == ReturnCode.Success
        @test 0 < sol.stats.njacs <= 20
    end
end

@testset "out-of-place long-horizon Jacobian reuse" begin
    f(u, p, t) = (du = similar(u); long_rober!(du, u, p, t); du)
    prob = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1.0e11))
    for alg in (TRBDF2(), KenCarp4())
        sol = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-10)
        @test sol.retcode == ReturnCode.Success
        @test 0 < sol.stats.njacs < 100
        @test sol.stats.njacs < sol.stats.nw
    end
end
