using OrdinaryDiffEqFIRK, LinearSolve, LinearAlgebra, SciMLOperators, Test
using SciMLOperators: AbstractSciMLOperator
using ADTypes: AutoForwardDiff, AutoFiniteDiff
using OrdinaryDiffEqFIRK: ComplexWOperator, set_W_gamma!

# 1D Brusselator: the diffusion coupling gives GMRES a system worth iterating on
const N_KRYLOV = 20

function brusselator_1d!(du, u, p, t)
    α = p[1]
    n = length(u) ÷ 2
    for i in 1:n
        x = u[i]
        y = u[n + i]
        xm = u[i == 1 ? n : i - 1]
        xp = u[i == n ? 1 : i + 1]
        ym = u[i == 1 ? 2n : n + i - 1]
        yp = u[i == n ? n + 1 : n + i + 1]
        du[i] = 1 + x^2 * y - 4x + α * (xm - 2x + xp)
        du[n + i] = 3x - x^2 * y + α * (ym - 2y + yp)
    end
    return nothing
end

u0_krylov = vcat(
    [1 + 0.5sinpi(2i / N_KRYLOV) for i in 1:N_KRYLOV],
    [3 + 0.3cospi(2i / N_KRYLOV) for i in 1:N_KRYLOV]
)
mm_krylov = Matrix(Diagonal(vcat(fill(1.0, N_KRYLOV), fill(2.0, N_KRYLOV))))
prob_krylov = ODEProblem(brusselator_1d!, u0_krylov, (0.0, 1.0), (10.0,))
prob_krylov_mm = ODEProblem(
    ODEFunction(brusselator_1d!; mass_matrix = mm_krylov), u0_krylov, (0.0, 1.0), (10.0,)
)
ref_krylov = solve(prob_krylov, RadauIIA5(), reltol = 1.0e-14, abstol = 1.0e-14).u[end]
ref_krylov_mm = solve(prob_krylov_mm, RadauIIA5(), reltol = 1.0e-14, abstol = 1.0e-14).u[end]

firk_algs = (RadauIIA3, RadauIIA5, RadauIIA9, AdaptiveRadau)

@testset "ComplexWOperator agrees with its dense form" begin
    n = 7
    J = randn(n, n)
    x = randn(ComplexF64, n)
    γ = 0.3 - 1.7im
    for mass_matrix in (I, 2.5I, Matrix(Diagonal(rand(n) .+ 1)), randn(n, n))
        W = ComplexWOperator(mass_matrix, complex(1.0), J, zeros(n))
        set_W_gamma!(W, γ)
        dense = -γ * (mass_matrix isa UniformScaling ? Matrix(mass_matrix, n, n) : mass_matrix) + J
        @test size(W) == (n, n)
        @test eltype(W) == ComplexF64
        @test W * x ≈ dense * x
        @test mul!(similar(x), W, x) ≈ dense * x
    end
end

@testset "$alg with a matrix-free linsolve" for alg in firk_algs
    integ = init(prob_krylov, alg(linsolve = KrylovJL_GMRES()))
    @test integ.cache.J isa AbstractSciMLOperator
    # RadauIIA3's single stage pair makes W1 itself the complex operator
    complex_W = if alg === RadauIIA3
        integ.cache.W1
    elseif alg === AdaptiveRadau
        integ.cache.W2[1]
    else
        integ.cache.W2
    end
    @test complex_W isa ComplexWOperator

    kry = solve(prob_krylov, alg(linsolve = KrylovJL_GMRES()), reltol = 1.0e-10, abstol = 1.0e-10)
    @test SciMLBase.successful_retcode(kry)
    # measured worst case over these methods is 4.8e-9; the margin is for CI variation
    @test norm(kry.u[end] - ref_krylov, Inf) < 1.0e-6

    mm = solve(
        prob_krylov_mm, alg(linsolve = KrylovJL_GMRES()), reltol = 1.0e-10, abstol = 1.0e-10
    )
    @test SciMLBase.successful_retcode(mm)
    @test norm(mm.u[end] - ref_krylov_mm, Inf) < 1.0e-5
end

@testset "$alg with concrete_jac" for alg in firk_algs
    # concrete J behind a Krylov solve: products go through the JVP, J is kept for preconditioning
    kry = solve(
        prob_krylov, alg(linsolve = KrylovJL_GMRES(), concrete_jac = true),
        reltol = 1.0e-10, abstol = 1.0e-10
    )
    @test SciMLBase.successful_retcode(kry)
    @test norm(kry.u[end] - ref_krylov, Inf) < 1.0e-6

    # concrete J behind a factorization: FIRK assembles and factorizes its own W
    lu = solve(prob_krylov, alg(concrete_jac = true), reltol = 1.0e-10, abstol = 1.0e-10)
    @test SciMLBase.successful_retcode(lu)
    @test norm(lu.u[end] - ref_krylov, Inf) < 1.0e-6
end

@testset "$alg with an operator jac_prototype" for alg in firk_algs
    n = 12
    A = -50 * Matrix(SymTridiagonal(fill(2.0, n), fill(-1.0, n - 1)))
    linear!(du, u, p, t) = mul!(du, A, u)
    prob = ODEProblem(
        ODEFunction(linear!; jac_prototype = MatrixOperator(copy(A))), ones(n), (0.0, 0.5)
    )
    sol = solve(prob, alg(linsolve = KrylovJL_GMRES()), reltol = 1.0e-10, abstol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)
    @test norm(sol.u[end] - exp(A * 0.5) * ones(n), Inf) < 1.0e-8

    @test_throws "requires a concrete matrix" solve(prob, alg())
end

@testset "GaussLegendre reports that it has no matrix-free path" begin
    @test_throws "does not support a matrix-free Jacobian" solve(
        prob_krylov, GaussLegendre(num_stages = 3, linsolve = KrylovJL_GMRES())
    )
end

@testset "$alg counts matrix-free Jacobian-vector products" for alg in firk_algs
    nf = Ref(0)
    counted!(du, u, p, t) = (nf[] += 1; brusselator_1d!(du, u, p, t))
    prob = ODEProblem(counted!, u0_krylov, (0.0, 1.0), (10.0,))
    for ad in (AutoForwardDiff(), AutoFiniteDiff())
        nf[] = 0
        sol = solve(
            prob, alg(linsolve = KrylovJL_GMRES(), autodiff = ad),
            reltol = 1.0e-8, abstol = 1.0e-8
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.stats.nf == nf[]
    end
end

@testset "$alg with a matrix-shaped state" for alg in firk_algs
    A = [-30.0 1.0; 2.0 -40.0]
    function matstate!(du, u, p, t)
        mul!(du, A, u)
        du .-= 0.1 .* u .^ 2
        return nothing
    end
    prob = ODEProblem(matstate!, [1.0 2.0 3.0; 4.0 5.0 6.0], (0.0, 0.2))
    ref = solve(prob, RadauIIA5(), reltol = 1.0e-13, abstol = 1.0e-13).u[end]
    sol = solve(prob, alg(linsolve = KrylovJL_GMRES()), reltol = 1.0e-10, abstol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)
    @test maximum(abs, sol.u[end] - ref) < 1.0e-8
end
