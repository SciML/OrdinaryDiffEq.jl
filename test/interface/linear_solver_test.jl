using Test, OrdinaryDiffEq

using SparseArrays, LinearSolve
using LinearAlgebra, Random
N = 30
AA = sprand(MersenneTwister(12), N, N, 0.5)
mm = sprand(MersenneTwister(123), N, N, 0.5)
A = MatrixOperator(AA)
M = MatrixOperator(mm'mm)
u0 = ones(N)
prob = ODEProblem(ODEFunction(A; mass_matrix = M), u0, (0.0, 1.0))

for alg in [Rosenbrock23(), Rosenbrock23(linsolve = KLUFactorization())]
    sol = solve(prob, alg)
    @test sol.stats.njacs == 0
    @test sol.stats.nw == 1
end

## OOP

foop(u, p, t) = jac(u, p, t) * u
jac(u, p, t) = spdiagm(0 => p)
paramjac(u, p, t) = SparseArrays.spdiagm(0 => u)

n = 2
p = collect(1.0:n)
u0 = ones(n)
tspan = [0.0, 1]
odef = ODEFunction(foop; jac = jac, jac_prototype = jac(u0, p, 0.0), paramjac = paramjac)

function g_helper(p; alg = Rosenbrock23(linsolve = LUFactorization()))
    prob = ODEProblem(odef, u0, tspan, p)
    soln = Array(solve(prob, alg; u0 = prob.u0, p = prob.p, abstol = 1e-4, reltol = 1e-4))[
        :, end]
    return soln
end
function g(p; kwargs...)
    soln = g_helper(p; kwargs...)
    return sum(soln)
end

@test isapprox(exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = KLUFactorization()));
    atol = 1e-3, rtol = 1e-3)
@test isapprox(exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = KrylovJL_GMRES()));
    atol = 1e-1, rtol = 1e-1)

@test isapprox(exp.(p), g_helper(p; alg = Rodas4(linsolve = KLUFactorization()));
    atol = 1e-3, rtol = 1e-3)
@test isapprox(exp.(p), g_helper(p; alg = Rodas4(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas4(linsolve = KrylovJL_GMRES())); atol = 1e-1,
    rtol = 1e-1)

@test isapprox(exp.(p), g_helper(p; alg = Rodas5(linsolve = KLUFactorization()));
    atol = 1e-3, rtol = 1e-3)
@test isapprox(exp.(p), g_helper(p; alg = Rodas5(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas5(linsolve = KrylovJL_GMRES())); atol = 1e-1,
    rtol = 1e-1)

@test isapprox(exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = KLUFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(
    exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = KrylovJL_GMRES()));
    atol = 1e-1, rtol = 1e-1)

@test isapprox(exp.(p), g_helper(p; alg = TRBDF2(linsolve = KLUFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(exp.(p), g_helper(p; alg = TRBDF2(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(
    exp.(p), g_helper(p; alg = TRBDF2(linsolve = KrylovJL_GMRES())); atol = 1e-1,
    rtol = 1e-1)

@test isapprox(exp.(p), g_helper(p; alg = KenCarp47(linsolve = KLUFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(exp.(p), g_helper(p; alg = KenCarp47(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(exp.(p), g_helper(p; alg = KenCarp47(linsolve = KrylovJL_GMRES()));
    atol = 1e-1, rtol = 1e-1)

@test isapprox(exp.(p), g_helper(p; alg = FBDF(linsolve = KrylovJL_GMRES()));
    atol = 1e-1, rtol = 1e-1)

## IIP

fiip(du, u, p, t) = du .= jac(u, p, t) * u
jac(du, u, p, t) = du .= spdiagm(0 => p)
paramjac(du, u, p, t) = du .= SparseArrays.spdiagm(0 => u)

n = 2
p = collect(1.0:n)
u0 = ones(n)
tspan = [0.0, 1]
odef = ODEFunction{true}(fiip; jac = jac, jac_prototype = jac(u0, p, 0.0),
    paramjac = paramjac)

function g_helper(p; alg = Rosenbrock23(linsolve = LUFactorization()))
    prob = ODEProblem(odef, u0, tspan, p)
    soln = Array(solve(prob, alg; u0 = prob.u0, p = prob.p, abstol = 1e-4, reltol = 1e-4))[
        :, end]
    return soln
end
function g(p; kwargs...)
    soln = g_helper(p; kwargs...)
    return sum(soln)
end

@test isapprox(exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = KLUFactorization()));
    atol = 1e-3, rtol = 1e-3)
@test isapprox(exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = KrylovJL_GMRES()));
    atol = 1e-1, rtol = 1e-1)

@test isapprox(exp.(p), g_helper(p; alg = Rodas4(linsolve = KLUFactorization()));
    atol = 1e-3, rtol = 1e-3)
@test isapprox(exp.(p), g_helper(p; alg = Rodas4(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas4(linsolve = KrylovJL_GMRES())); atol = 1e-1,
    rtol = 1e-1)

@test isapprox(exp.(p), g_helper(p; alg = Rodas5(linsolve = KLUFactorization()));
    atol = 1e-3, rtol = 1e-3)
@test isapprox(exp.(p), g_helper(p; alg = Rodas5(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas5(linsolve = KrylovJL_GMRES())); atol = 1e-1,
    rtol = 1e-1)

@test isapprox(exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = KLUFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(
    exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = KrylovJL_GMRES()));
    atol = 1e-1, rtol = 1e-1)

@test isapprox(exp.(p), g_helper(p; alg = TRBDF2(linsolve = KLUFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(exp.(p), g_helper(p; alg = TRBDF2(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(
    exp.(p), g_helper(p; alg = TRBDF2(linsolve = KrylovJL_GMRES())); atol = 1e-1,
    rtol = 1e-1)

@test isapprox(exp.(p), g_helper(p; alg = KenCarp47(linsolve = KLUFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(exp.(p), g_helper(p; alg = KenCarp47(linsolve = UMFPACKFactorization()));
    atol = 1e-1, rtol = 1e-1)
@test isapprox(exp.(p), g_helper(p; alg = KenCarp47(linsolve = KrylovJL_GMRES()));
    atol = 1e-1, rtol = 1e-1)

using OrdinaryDiffEq, StaticArrays, LinearSolve, ParameterizedFunctions

hires = @ode_def Hires begin
    dy1 = -1.71f0 * y1 + 0.43f0 * y2 + 8.32f0 * y3 + 0.0007f0 + 1.0f-18 * t
    dy2 = 1.71f0 * y1 - 8.75f0 * y2
    dy3 = -10.03f0 * y3 + 0.43f0 * y4 + 0.035f0 * y5
    dy4 = 8.32f0 * y2 + 1.71f0 * y3 - 1.12f0 * y4
    dy5 = -1.745f0 * y5 + 0.43f0 * y6 + 0.43f0 * y7
    dy6 = -280.0f0 * y6 * y8 + 0.69f0 * y4 + 1.71f0 * y5 - 0.43f0 * y6 + 0.69f0 * y7
    dy7 = 280.0f0 * y6 * y8 - 1.81f0 * y7
    dy8 = -280.0f0 * y6 * y8 + 1.81f0 * y7
end

u0 = zeros(8)
u0[1] = 1
u0[8] = 0.0057

probiip = ODEProblem{true}(hires, u0, (0.0, 10.0))
proboop = ODEProblem{false}(hires, u0, (0.0, 10.0))
probstatic = ODEProblem{false}(hires, SVector{8}(u0), (0.0, 10.0))
probiipf32 = ODEProblem{true}(hires, Float32.(u0), (0.0f0, 10.0f0))
proboopf32 = ODEProblem{false}(hires, Float32.(u0), (0.0f0, 10.0f0))
probstaticf32 = ODEProblem{false}(hires, SVector{8}(Float32.(u0)), (0.0f0, 10.0f0))
probs = (; probiip, proboop, probstatic)
probsf32 = (; probiipf32, proboopf32, probstaticf32)
qndf = QNDF()
krylov_qndf = QNDF(linsolve = KrylovJL_GMRES())
fbdf = FBDF()
krylov_fbdf = FBDF(linsolve = KrylovJL_GMRES())
rodas = Rodas5P()
krylov_rodas = Rodas5P(linsolve = KrylovJL_GMRES())
solvers = (; qndf, krylov_qndf, rodas, krylov_rodas, fbdf, krylov_fbdf)

refsol = solve(probiip, FBDF(), abstol = 1e-12, reltol = 1e-12)
@testset "Hires calc_W tests" begin
    @testset "$probname" for (probname, prob) in pairs(probs)
        @testset "$solname" for (solname, solver) in pairs(solvers)
            sol = solve(prob, solver, abstol = 1e-12, reltol = 1e-12, maxiters = 2e4)
            @test sol.retcode == ReturnCode.Success
            @test isapprox(sol.u[end], refsol.u[end], rtol = 2e-8, atol = 1e-10)
        end
    end
end

@testset "Hires Float32 calc_W tests" begin
    @testset "$probname" for (probname, prob) in pairs(probsf32)
        @testset "$solname" for (solname, solver) in pairs(solvers)
            sol = solve(prob, solver, maxiters = 2e4)
            @test sol.retcode == ReturnCode.Success
            @test isapprox(sol.u[end], refsol.u[end], rtol = 5e-3, atol = 1e-6)
        end
    end
end
