using OrdinaryDiffEqDifferentiation: rhs_evals_per_jvp, base_evals_per_jvp_point, JVPCache
using OrdinaryDiffEqBDF, LinearSolve, LinearAlgebra, ADTypes, SciMLBase, Test
using SciMLOperators: MatrixOperator, update_coefficients!

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
    @test rhs_evals_per_jvp(AutoFiniteDiff()) == 1
    @test rhs_evals_per_jvp(AutoFiniteDiff(fdjtype = Val(:central))) == 1
    @test rhs_evals_per_jvp(AutoSparse(AutoFiniteDiff())) == 1
    @test rhs_evals_per_jvp(AutoSparse(AutoForwardDiff())) == 1

    @test base_evals_per_jvp_point(AutoFiniteDiff()) == 1
    @test base_evals_per_jvp_point(AutoForwardDiff()) == 0
end

# The base evaluation `f(u)` is fixed by `update_coefficients!` and reused by every product
# taken at that point, so a run of `K` products costs `K + 1` evaluations rather than `2K`.
@testset "the finite-difference base evaluation is shared across products" begin
    K = 5
    prob = matrix_free
    @testset "$adname" for (adname, ad, per_point) in (
            ("AutoFiniteDiff", AutoFiniteDiff(), 1), ("AutoForwardDiff", AutoForwardDiff(), 0),
        )
        J = JVPCache(prob.f, similar(u0), u0, prob.p, 0.0; autodiff = ad)
        v, Jv = ones(N), zeros(N)

        NCALLS[] = 0
        update_coefficients!(J, u0, prob.p, 0.0)
        for _ in 1:K
            mul!(Jv, J, v)
        end
        @test NCALLS[] == K + per_point
        @test J.njvps == K
        @test J.nbase_evals == per_point

        # A second point costs its own base evaluation, and only one however many times it
        # is set before a product is taken.
        NCALLS[] = 0
        update_coefficients!(J, u0 .+ 1, prob.p, 0.0)
        update_coefficients!(J, u0 .+ 1, prob.p, 0.0)
        mul!(Jv, J, v)
        @test NCALLS[] == 1 + per_point

        # The JVP is still the directional derivative it was before any of this caching.
        fd = (similar(u0), similar(u0))
        h = 1.0e-6
        allencahn!(fd[1], u0 .+ 1 .+ h .* v, prob.p, 0.0)
        allencahn!(fd[2], u0 .+ 1, prob.p, 0.0)
        @test Jv ≈ (fd[1] .- fd[2]) ./ h rtol = 1.0e-4
    end
end
