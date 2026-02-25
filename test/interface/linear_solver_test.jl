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
    soln = Array(solve(prob, alg; u0 = prob.u0, p = prob.p, abstol = 1.0e-4, reltol = 1.0e-4))[
        :, end,
    ]
    return soln
end
function g(p; kwargs...)
    soln = g_helper(p; kwargs...)
    return sum(soln)
end

@test isapprox(
    exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = KLUFactorization()));
    atol = 1.0e-2, rtol = 1.0e-2
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = KrylovJL_GMRES()));
    atol = 1.0e-1, rtol = 1.0e-1
)

@test isapprox(
    exp.(p), g_helper(p; alg = Rodas4(linsolve = KLUFactorization()));
    atol = 1.0e-3, rtol = 1.0e-3
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas4(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas4(linsolve = KrylovJL_GMRES())); atol = 1.0e-1,
    rtol = 1.0e-1
)

@test isapprox(
    exp.(p), g_helper(p; alg = Rodas5(linsolve = KLUFactorization()));
    atol = 1.0e-3, rtol = 1.0e-3
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas5(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas5(linsolve = KrylovJL_GMRES())); atol = 1.0e-1,
    rtol = 1.0e-1
)

@test isapprox(
    exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = KLUFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = KrylovJL_GMRES()));
    atol = 1.0e-1, rtol = 1.0e-1
)

@test isapprox(
    exp.(p), g_helper(p; alg = TRBDF2(linsolve = KLUFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = TRBDF2(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = TRBDF2(linsolve = KrylovJL_GMRES())); atol = 1.0e-1,
    rtol = 1.0e-1
)

@test isapprox(
    exp.(p), g_helper(p; alg = KenCarp47(linsolve = KLUFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = KenCarp47(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = KenCarp47(linsolve = KrylovJL_GMRES()));
    atol = 1.0e-1, rtol = 1.0e-1
)

@test isapprox(
    exp.(p), g_helper(p; alg = FBDF(linsolve = KrylovJL_GMRES()));
    atol = 1.0e-1, rtol = 1.0e-1
)

## IIP

fiip(du, u, p, t) = du .= jac(u, p, t) * u
jac(du, u, p, t) = du .= spdiagm(0 => p)
paramjac(du, u, p, t) = du .= SparseArrays.spdiagm(0 => u)

n = 2
p = collect(1.0:n)
u0 = ones(n)
tspan = [0.0, 1]
odef = ODEFunction{true}(
    fiip; jac = jac, jac_prototype = jac(u0, p, 0.0),
    paramjac = paramjac
)

function g_helper(p; alg = Rosenbrock23(linsolve = LUFactorization()))
    prob = ODEProblem(odef, u0, tspan, p)
    soln = Array(solve(prob, alg; u0 = prob.u0, p = prob.p, abstol = 1.0e-4, reltol = 1.0e-4))[
        :, end,
    ]
    return soln
end
function g(p; kwargs...)
    soln = g_helper(p; kwargs...)
    return sum(soln)
end

@test isapprox(
    exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = KLUFactorization()));
    atol = 1.0e-2, rtol = 1.0e-2
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rosenbrock23(linsolve = KrylovJL_GMRES()));
    atol = 1.0e-1, rtol = 1.0e-1
)

@test isapprox(
    exp.(p), g_helper(p; alg = Rodas4(linsolve = KLUFactorization()));
    atol = 1.0e-3, rtol = 1.0e-3
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas4(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas4(linsolve = KrylovJL_GMRES())); atol = 1.0e-1,
    rtol = 1.0e-1
)

@test isapprox(
    exp.(p), g_helper(p; alg = Rodas5(linsolve = KLUFactorization()));
    atol = 1.0e-3, rtol = 1.0e-3
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas5(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = Rodas5(linsolve = KrylovJL_GMRES())); atol = 1.0e-1,
    rtol = 1.0e-1
)

@test isapprox(
    exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = KLUFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = ImplicitEuler(linsolve = KrylovJL_GMRES()));
    atol = 1.0e-1, rtol = 1.0e-1
)

@test isapprox(
    exp.(p), g_helper(p; alg = TRBDF2(linsolve = KLUFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = TRBDF2(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = TRBDF2(linsolve = KrylovJL_GMRES())); atol = 1.0e-1,
    rtol = 1.0e-1
)

@test isapprox(
    exp.(p), g_helper(p; alg = KenCarp47(linsolve = KLUFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = KenCarp47(linsolve = UMFPACKFactorization()));
    atol = 1.0e-1, rtol = 1.0e-1
)
@test isapprox(
    exp.(p), g_helper(p; alg = KenCarp47(linsolve = KrylovJL_GMRES()));
    atol = 1.0e-1, rtol = 1.0e-1
)

using OrdinaryDiffEq, StaticArrays, LinearSolve

# In-place version
function hires!(du, u, p, t)
    T = eltype(u)
    du[1] = T(-1.71) * u[1] + T(0.43) * u[2] + T(8.32) * u[3] + T(0.0007) + T(1.0e-18) * t
    du[2] = T(1.71) * u[1] - T(8.75) * u[2]
    du[3] = T(-10.03) * u[3] + T(0.43) * u[4] + T(0.035) * u[5]
    du[4] = T(8.32) * u[2] + T(1.71) * u[3] - T(1.12) * u[4]
    du[5] = T(-1.745) * u[5] + T(0.43) * u[6] + T(0.43) * u[7]
    du[6] = T(-280.0) * u[6] * u[8] + T(0.69) * u[4] + T(1.71) * u[5] - T(0.43) * u[6] +
        T(0.69) * u[7]
    du[7] = T(280.0) * u[6] * u[8] - T(1.81) * u[7]
    return du[8] = T(-280.0) * u[6] * u[8] + T(1.81) * u[7]
end

# Out-of-place version
function hires(u, p, t)
    T = eltype(u)
    du1 = T(-1.71) * u[1] + T(0.43) * u[2] + T(8.32) * u[3] + T(0.0007) + T(1.0e-18) * t
    du2 = T(1.71) * u[1] - T(8.75) * u[2]
    du3 = T(-10.03) * u[3] + T(0.43) * u[4] + T(0.035) * u[5]
    du4 = T(8.32) * u[2] + T(1.71) * u[3] - T(1.12) * u[4]
    du5 = T(-1.745) * u[5] + T(0.43) * u[6] + T(0.43) * u[7]
    du6 = T(-280.0) * u[6] * u[8] + T(0.69) * u[4] + T(1.71) * u[5] - T(0.43) * u[6] +
        T(0.69) * u[7]
    du7 = T(280.0) * u[6] * u[8] - T(1.81) * u[7]
    du8 = T(-280.0) * u[6] * u[8] + T(1.81) * u[7]

    if u isa SVector
        return SVector(du1, du2, du3, du4, du5, du6, du7, du8)
    else
        return [du1, du2, du3, du4, du5, du6, du7, du8]
    end
end

u0 = zeros(8)
u0[1] = 1
u0[8] = 0.0057

probiip = ODEProblem{true}(hires!, u0, (0.0, 10.0))
proboop = ODEProblem{false}(hires, u0, (0.0, 10.0))
probstatic = ODEProblem{false}(hires, SVector{8}(u0), (0.0, 10.0))
probiipf32 = ODEProblem{true}(hires!, Float32.(u0), (0.0f0, 10.0f0))
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

refsol = solve(probiip, FBDF(), abstol = 1.0e-12, reltol = 1.0e-12)
@testset "Hires calc_W tests" begin
    @testset "$probname" for (probname, prob) in pairs(probs)
        @testset "$solname" for (solname, solver) in pairs(solvers)
            sol = solve(prob, solver, abstol = 1.0e-12, reltol = 1.0e-12, maxiters = 2.0e4)
            @test sol.retcode == ReturnCode.Success
            @test isapprox(sol.u[end], refsol.u[end], rtol = 1.0e-8, atol = 1.0e-10)
        end
    end
end

@testset "Hires Float32 calc_W tests" begin
    @testset "$probname" for (probname, prob) in pairs(probsf32)
        @testset "$solname" for (solname, solver) in pairs(solvers)
            sol = solve(prob, solver, maxiters = 2.0e4)
            @test sol.retcode == ReturnCode.Success
            @test isapprox(sol.u[end], refsol.u[end], rtol = 2.0e-3, atol = 1.0e-6)
        end
    end
end
