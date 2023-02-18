# Solve the Landau-Zener problem i ψ' = H(t) ψ, with H(t) = [t 1;1 -t]

using Test
using StaticArrays, LinearAlgebra
using OrdinaryDiffEq, DiffEqBase

H(t) = -im * (@SMatrix [t 1; 1 -t])

fun(ψ, p, t) = oftype(ψ, H(t) * ψ)
fun_inplace(dψ, ψ, p, t) = (dψ .= H(t) * ψ)

T = 0.1
tspan = (0, T)
explicit = [Midpoint, RK4, DP5, Tsit5, Vern7]
implicit = [ImplicitEuler, Trapezoid, Kvaerno3, Rosenbrock23]

@testset "Complex Tests on Explicit Methods. alg=$alg" for alg in explicit
    for f in (fun, fun_inplace)
        ψ0 = [1.0 + 0.0im; 0.0]
        prob = ODEProblem(f, ψ0, (-T, T))
        sol = solve(prob, alg())
        @test norm(sol(T))≈1 atol=1e-2
    end
    ψ0 = @SArray [1.0 + 0.0im; 0.0]
    prob = ODEProblem(fun, ψ0, (-T, T))
    sol = solve(prob, alg())
    @test norm(sol(T))≈1 atol=1e-2
end

@testset "Complex Tests on Implicit Autodiff Methods. alg=$alg" for alg in implicit
    @test_broken begin
        for f in (fun, fun_inplace)
            ψ0 = [1.0 + 0.0im; 0.0]
            prob = ODEProblem(f, ψ0, (-T, T))
            sol = solve(prob, alg())
            @test norm(sol(T))≈1 atol=1e-2
        end
        ψ0 = @SArray [1.0 + 0.0im; 0.0]
        prob = ODEProblem(fun, ψ0, (-T, T))
        sol = solve(prob, alg())
        @test norm(sol(T))≈1 atol=1e-2
    end
end

@testset "Complex Tests on Implicit Finite Diff Methods. alg=$alg" for alg in implicit
    ψ0 = [1.0 + 0.0im; 0.0]
    prob = ODEProblem(fun_inplace, ψ0, (-T, T))
    sol = solve(prob, alg(autodiff = false))
    @test norm(sol(T))≈1 atol=1e-2
end

@testset "Complex Tests on Implicit Finite Diff Out-of-place Methods. alg=$alg" for alg in implicit
    ψ0 = [1.0 + 0.0im; 0.0]
    prob = ODEProblem(fun, ψ0, (-T, T))
    sol = solve(prob, alg(autodiff = false))
    @test norm(sol(T))≈1 atol=1e-2
end

@testset "Complex Tests on Implicit Finite Diff Out-of-place Methods SArray. alg=$alg" for alg in implicit
    ψ0 = @SArray [1.0 + 0.0im; 0.0]
    prob = ODEProblem(fun, ψ0, (-T, T))
    sol = solve(prob, alg(autodiff = false))
    @test norm(sol(T))≈1 atol=1e-2
end
