using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using NonlinearSolve: NewtonRaphson
using LinearSolve, LinearAlgebra, SparseArrays, ADTypes
using SciMLBase
using SciMLOperators: MatrixOperator
using Test

# 1D Brusselator with a hand-assembled sparse Jacobian pattern (tridiagonal blocks).
const N_b = 8
function bruss1d!(du, u, p, t)
    A, B, alpha, dx = p
    a = alpha / dx^2
    n = N_b
    @inbounds for i in 1:n
        im1 = i == 1 ? n : i - 1
        ip1 = i == n ? 1 : i + 1
        x = u[i]
        y = u[n + i]
        du[i] = a * (u[im1] + u[ip1] - 2x) + B + x^2 * y - (A + 1) * x
        du[n + i] = a * (u[n + im1] + u[n + ip1] - 2y) + A * x - x^2 * y
    end
    return nothing
end
function bruss1d_jacproto(n)
    Is = Int[]
    Js = Int[]
    for i in 1:n
        im1 = i == 1 ? n : i - 1
        ip1 = i == n ? 1 : i + 1
        for (r, c) in (
                (i, im1), (i, i), (i, ip1), (i, n + i),
                (n + i, n + im1), (n + i, n + i), (n + i, n + ip1), (n + i, i),
            )
            push!(Is, r)
            push!(Js, c)
        end
    end
    return sparse(Is, Js, ones(length(Is)), 2n, 2n)
end
u0 = vcat(
    [22 * (i / (N_b + 1) * (1 - i / (N_b + 1)))^(3 / 2) for i in 1:N_b],
    [27 * (i / (N_b + 1) * (1 - i / (N_b + 1)))^(3 / 2) for i in 1:N_b]
)
p = (3.4, 1.0, 10.0, 1 / N_b)
f_sparse = ODEFunction(bruss1d!; jac_prototype = bruss1d_jacproto(N_b))
prob = ODEProblem(f_sparse, u0, (0.0, 2.0), p)

refsol = solve(prob, FBDF(); reltol = 1.0e-10, abstol = 1.0e-12)

nsa = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff()))

@testset "NSA W-reuse keeps inner Jacobian sparse" begin
    integ = init(prob, FBDF(nlsolve = nsa); reltol = 1.0e-8, abstol = 1.0e-10)
    nsacache = integ.cache.nlsolver.cache
    @test nsacache.W isa SparseMatrixCSC
    # The inner NonlinearSolve Jacobian is the reused W as an operator; its underlying
    # matrix must mirror W's structure (stay sparse), not densify.
    J = integ.cache.nlsolver.cache.cache.jac_cache.J
    @test J isa MatrixOperator
    @test convert(AbstractMatrix, J) isa SparseMatrixCSC
end

@testset "NSA + sparse-only linsolve (KLU) solves" begin
    for alg in (
            FBDF(linsolve = KLUFactorization(), nlsolve = nsa),
            TRBDF2(linsolve = KLUFactorization(), nlsolve = nsa),
        )
        sol = solve(prob, alg; reltol = 1.0e-8, abstol = 1.0e-10)
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end] ≈ refsol.u[end] rtol = 1.0e-5
    end
end

# Resizing the integrator must rebuild the inner Jacobian at the new size: the jac
# closure and jac_prototype are fixed at build time, so the length-mismatch path in
# initialize! has to refresh both against the resized W (this crashed with
# `DimensionMismatch: B has leading dimension 2, but needs 1` when jac_prototype
# followed the W structure).
@testset "resize with W-reuse rebuilds inner Jacobian" begin
    fgrow = function (du, u, p, t)
        for i in 1:length(u)
            du[i] = (0.3 / length(u)) * u[i]
        end
        return nothing
    end
    condition = (u, t, integrator) -> 1 - maximum(u)
    affect! = function (integrator)
        u = integrator.u
        maxidx = findmax(u)[2]
        resize!(integrator, length(u) + 1)
        Θ = 0.3
        u[maxidx] = Θ
        u[end] = 1 - Θ
        return nothing
    end
    cb = ContinuousCallback(condition, affect!)
    probg = ODEProblem(fgrow, [0.2], (0.0, 10.0))
    # NOTE: assert the resize worked and nothing crashed (mirroring
    # test/Integrators_II/ode_cache_tests.jl). Whether the marginally-stable
    # grow-a-cell dynamics finishes with Success is method- and
    # rounding-sensitive (ImplicitEuler+NLNewton is Unstable on it too), so
    # retcode is only checked for TRBDF2.
    for alg in (
            TRBDF2(nlsolve = nsa), KenCarp4(nlsolve = nsa),
            ImplicitEuler(nlsolve = nsa),
        )
        sol = solve(probg, alg; callback = cb, dt = 1 / 2)
        @test length(sol.u[end]) > 1
    end
    sol = solve(probg, TRBDF2(nlsolve = nsa); callback = cb, dt = 1 / 2)
    @test SciMLBase.successful_retcode(sol)
end
