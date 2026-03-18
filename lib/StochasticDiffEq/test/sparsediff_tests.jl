using Test
using StochasticDiffEq
using SparseArrays
using LinearAlgebra

function f(dx, x, p, t)
    for i in 2:(length(x) - 1)
        dx[i] = x[i - 1] - 2x[i] + x[i + 1]
    end
    dx[1] = -2x[1] + x[2]
    dx[end] = x[end - 1] - 2x[end]
    return nothing
end

function g(dx, x, p, t)
    dx .= 0
    return nothing
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
colorvec = repeat(1:3, 10)[1:10]
u0 = [1.0, 2.0, 3, 4, 5, 5, 4, 3, 2, 1]
tspan = (0.0, 10.0)
sdefun_sp = SDEFunction(f, g, colorvec = colorvec, jac_prototype = jac_sp)
prob_sp = SDEProblem(sdefun_sp, u0, tspan)
prob_std = SDEProblem(f, g, u0, tspan)

sol_sp = solve(prob_sp, SKenCarp(autodiff = false))
@test sol_sp.retcode == ReturnCode.Success #test sparse finitediff
sol = solve(prob_std, SKenCarp(autodiff = false))
@test sol_sp.u[end] ≈ sol.u[end] atol = 1.0e-4
@test length(sol_sp.t) == length(sol.t)

sol_sp = solve(prob_sp, SKenCarp())
sol = solve(prob_std, SKenCarp())
@test sol_sp.u[end] ≈ sol.u[end] atol = 1.0e-4
@test length(sol_sp.t) == length(sol.t)

#sol_sp=solve(prob_sp,SKenCarp(autodiff=false),abstol=1e-10,reltol=1e-10)
#sol=solve(prob_std,SKenCarp(autodiff=false),abstol=1e-10,reltol=1e-10)
#@test sol_sp.u[end]≈sol.u[end]
#@test length(sol_sp.t)==length(sol.t)
