using OrdinaryDiffEqDifferentiation: rhs_evals_per_jvp
using OrdinaryDiffEqBDF, LinearSolve, LinearAlgebra, ADTypes, SciMLBase, Test
using SciMLOperators: MatrixOperator

# `stats.nf` counts calls to `f`. Every configuration below counts them for real and
# checks the solver reported exactly that many.
const N = 24
const dx = 1 / (N + 1)
const NCALLS = Ref(0)
function allencahn!(du, u, p, t)
    NCALLS[] += 1
    for i in 1:N
        um = i == 1 ? -one(eltype(u)) : u[i - 1]
        up = i == N ? one(eltype(u)) : u[i + 1]
        du[i] = 1.0e-3 * (um - 2u[i] + up) / dx^2 + u[i] - u[i]^3
    end
    return nothing
end
function allencahn_jac!(J, u, p, t)
    fill!(J, 0)
    for i in 1:N
        J[i, i] = -2.0e-3 / dx^2 + 1 - 3u[i]^2
        i > 1 && (J[i, i - 1] = 1.0e-3 / dx^2)
        i < N && (J[i, i + 1] = 1.0e-3 / dx^2)
    end
    return J
end
u0 = [tanh((i * dx - 0.5) / 0.1) for i in 1:N]
tspan = (0.0, 1.0)

matrix_free = ODEProblem(allencahn!, u0, tspan)
operator_jac = ODEProblem(
    ODEFunction(
        allencahn!; jac = allencahn_jac!,
        jac_prototype = MatrixOperator(zeros(N, N); update_func! = allencahn_jac!)
    ), u0, tspan
)

function counted_nf(prob, alg)
    solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-6, save_everystep = false)
    NCALLS[] = 0
    sol = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-6, save_everystep = false)
    @test SciMLBase.successful_retcode(sol)
    return NCALLS[], sol.stats.nf
end

@testset "stats.nf counts every RHS evaluation the Krylov solve costs" begin
    @testset "$adname" for (adname, ad) in
        (("AutoForwardDiff", AutoForwardDiff()), ("AutoFiniteDiff", AutoFiniteDiff()))
        # Fully matrix-free: the JVP operator performs every product.
        real_f, nf = counted_nf(matrix_free, FBDF(linsolve = KrylovJL_GMRES(), autodiff = ad))
        @test nf == real_f

        # Concrete J kept for the preconditioner, products still done by the JVP operator.
        real_f, nf = counted_nf(
            matrix_free, FBDF(linsolve = KrylovJL_GMRES(), concrete_jac = true, autodiff = ad)
        )
        @test nf == real_f

        # A user-supplied `MatrixOperator` Jacobian applies a stored matrix: no `f` calls.
        real_f, nf = counted_nf(operator_jac, FBDF(linsolve = KrylovJL_GMRES(), autodiff = ad))
        @test nf == real_f

        # Direct solver, no products at all.
        real_f, nf = counted_nf(matrix_free, FBDF(autodiff = ad))
        @test nf == real_f
    end
end

@testset "rhs_evals_per_jvp" begin
    @test rhs_evals_per_jvp(AutoForwardDiff()) == 1
    @test rhs_evals_per_jvp(AutoFiniteDiff()) == 2
    @test rhs_evals_per_jvp(AutoFiniteDiff(fdjtype = Val(:central))) == 2
    @test rhs_evals_per_jvp(AutoSparse(AutoFiniteDiff())) == 2
    @test rhs_evals_per_jvp(AutoSparse(AutoForwardDiff())) == 1
end
