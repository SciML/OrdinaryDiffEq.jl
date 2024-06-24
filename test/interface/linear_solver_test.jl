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
