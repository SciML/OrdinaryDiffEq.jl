using Test
using OrdinaryDiffEq
using SparseArrays
using LinearAlgebra
using LinearSolve
import DifferentiationInterface as DI
using SparseConnectivityTracer
using SparseMatrixColorings
using ADTypes

if isempty(VERSION.prerelease)
    using Enzyme
end

## in-place
#https://github.com/JuliaDiffEq/SparseDiffTools.jl/blob/master/test/test_integration.jl
function f_ip(dx, x, p, t)
    for i in 2:(length(x) - 1)
        dx[i] = x[i - 1] - 2x[i] + x[i + 1]
    end
    dx[1] = -2x[1] + x[2]
    dx[end] = x[end - 1] - 2x[end]
    return nothing
end

## out-of-place
function f_oop(x, p, t)
    dx = similar(x)
    for i in 2:(length(x) - 1)
        dx[i] = x[i - 1] - 2x[i] + x[i + 1]
    end
    dx[1] = -2x[1] + x[2]
    dx[end] = x[end - 1] - 2x[end]
    return dx
end

function second_derivative_stencil(N)
    A = zeros(N, N)
    for i in 1:N, j in 1:N
        (j - i == -1 || j - i == 1) && (A[i, j] = 1)
        j - i == 0 && (A[i, j] = -2)
    end
    return A
end

function generate_sparsity_pattern(N::Integer)
    dl = repeat([1.0], N - 1)
    du = repeat([1.0], N - 1)
    d = repeat([-2.0], N)
    return Tridiagonal(dl, d, du)
end

jac_sp = sparse(generate_sparsity_pattern(10))
#jac = second_derivative_stencil(10)
colors = repeat(1:3, 10)[1:10]
u0 = [1.0, 2.0, 3, 4, 5, 5, 4, 3, 2, 1]
tspan = (0.0, 10.0)

adchoices = if isempty(VERSION.prerelease)
    [
        AutoForwardDiff(), AutoFiniteDiff(),
        AutoEnzyme(mode = Enzyme.Forward, function_annotation = Enzyme.Const),
    ]
else
    [AutoForwardDiff(), AutoFiniteDiff()]
end

for f in [f_oop, f_ip]
    odefun_std = ODEFunction(f)
    prob_std = ODEProblem(odefun_std, u0, tspan)

    for ad in adchoices, linsolve in [nothing, LinearSolve.KrylovJL_GMRES()]
        for Solver in [Rodas5, Rosenbrock23, Trapezoid, KenCarp4, FBDF]
            for tol in [nothing, 1.0e-10]
                sol_std = solve(prob_std, Solver(autodiff = ad, linsolve = linsolve), reltol = tol, abstol = tol)
                @test sol_std.retcode == ReturnCode.Success
                for (i, prob) in enumerate(
                        map(
                            f -> ODEProblem(f, u0, tspan),
                            [
                                ODEFunction(
                                    f, colorvec = colors,
                                    jac_prototype = jac_sp
                                ),
                                ODEFunction(
                                    f, colorvec = colors,
                                    jac_prototype = jac_sp, mass_matrix = I(length(u0))
                                ),
                                ODEFunction(f, jac_prototype = jac_sp),
                                ODEFunction(
                                    f, colorvec = colors,
                                    sparsity = jac_sp
                                ),
                            ]
                        )
                    )
                    sol = solve(prob, Solver(autodiff = ad, linsolve = linsolve), reltol = tol, abstol = tol)
                    @test sol.retcode == ReturnCode.Success
                    if tol != nothing
                        @test sol_std.u[end] ≈ sol.u[end] atol = tol
                    else
                        @test sol_std.u[end] ≈ sol.u[end]
                    end
                    @test length(sol_std.t) == length(sol.t)
                end
            end
        end
    end
end

# test for https://github.com/SciML/OrdinaryDiffEq.jl/issues/2653#issuecomment-2778430025

using LinearAlgebra, SparseArrays
using OrdinaryDiffEq

function f(du, u, p, t)
    du[1] = u[1]
    return du
end

function jac(J::SparseMatrixCSC, u, p, t)
    @assert nnz(J) == 1  # mirrors the strict behavior of SparseMatrixColorings
    nonzeros(J)[1] = 1
    return J
end

u0 = ones(10)
jac_prototype = sparse(Diagonal(vcat(1, zeros(9))))

fun = ODEFunction(f; jac, jac_prototype)
prob = ODEProblem(fun, u0, (0.0, 1.0))
@test_nowarn sol = solve(prob, Rodas4(); reltol = 1.0e-8, abstol = 1.0e-8)
