using OrdinaryDiffEqDifferentiation: set_linear_reltol!
using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, OrdinaryDiffEqRosenbrock
using LinearSolve, LinearAlgebra, Test

# 1D Allen-Cahn: stiff enough that the Newton path takes several Krylov iterations.
const N = 24
const dx = 1 / (N + 1)
function allencahn!(du, u, p, t)
    for i in 1:N
        um = i == 1 ? -one(eltype(u)) : u[i - 1]
        up = i == N ? one(eltype(u)) : u[i + 1]
        du[i] = 1.0e-3 * (um - 2u[i] + up) / dx^2 + u[i] - u[i]^3
    end
    return nothing
end
u0 = [tanh((i * dx - 0.5) / 0.1) for i in 1:N]
prob = ODEProblem(allencahn!, u0, (0.0, 0.5))

newton_cache(integrator) = integrator.cache.nlsolver.cache.linsolve

@testset "integrator tolerance reaches the Krylov cache" begin
    # Newton path: `build_nlsolver` builds the `LinearCache` without tolerances, so
    # `dolinsolve` is the only thing that can put the integrator's `reltol` there.
    # A Newton direction is floored at `sqrt(eps)` however loose the integrator is, because
    # a correction solved only to `reltol` cannot drive the stage to `reltol` (#4293), so a
    # looser tolerance clamps and a tighter one passes through unchanged.
    for tol in (1.0e-6, 1.0e-10)
        integrator = init(prob, FBDF(linsolve = KrylovJL_GMRES()); reltol = tol, abstol = tol)
        for _ in 1:5
            step!(integrator)
        end
        @test newton_cache(integrator).reltol == min(tol, sqrt(eps(Float64)))
    end

    integrator = init(prob, KenCarp4(linsolve = KrylovJL_GMRES()); reltol = 1.0e-7, abstol = 1.0e-7)
    for _ in 1:5
        step!(integrator)
    end
    @test newton_cache(integrator).reltol == 1.0e-7

    # Rosenbrock already sets both tolerances at `init`; nothing may undo that.
    integrator = init(prob, Rodas5P(linsolve = KrylovJL_GMRES()); reltol = 1.0e-7, abstol = 1.0e-7)
    for _ in 1:5
        step!(integrator)
    end
    @test integrator.cache.linsolve.reltol == 1.0e-7
end

@testset "set_linear_reltol! only touches solvers that have a tolerance" begin
    A = Matrix(SymTridiagonal(fill(2.0, 8), fill(-1.0, 7)))
    b = ones(8)

    krylov = init(LinearProblem(copy(A), copy(b)), KrylovJL_GMRES())
    @test set_linear_reltol!(krylov, 1.0e-5).reltol == 1.0e-5

    # A factorization has no tolerance to update; `update_tolerances!` throws for one,
    # so the guard must keep it out.
    lu = init(LinearProblem(copy(A), copy(b)), LUFactorization())
    before = lu.reltol
    @test set_linear_reltol!(lu, 1.0e-5).reltol == before

    # Per-component tolerances have no scalar to hand a `LinearCache`.
    vec_tol = init(LinearProblem(copy(A), copy(b)), KrylovJL_GMRES())
    before = vec_tol.reltol
    @test set_linear_reltol!(vec_tol, fill(1.0e-5, 8)).reltol == before
end
