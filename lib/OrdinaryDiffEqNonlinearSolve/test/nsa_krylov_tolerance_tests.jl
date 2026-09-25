using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg, _inner_lintol
using NonlinearSolve, NonlinearSolveBase
using LinearSolve, LinearAlgebra, ADTypes, SciMLBase
using Test

# 1D Allen-Cahn: stiff enough that the Newton path takes several Krylov iterations, and
# small enough to build a dense reference from.
const N = 24
const DX = 1 / (N + 1)
function allencahn!(du, u, p, t)
    for i in 1:N
        um = i == 1 ? -one(eltype(u)) : u[i - 1]
        up = i == N ? one(eltype(u)) : u[i + 1]
        du[i] = 1.0e-3 * (um - 2u[i] + up) / DX^2 + u[i] - u[i]^3
    end
    return nothing
end
u0 = [tanh((i * DX - 0.5) / 0.1) for i in 1:N]
prob = ODEProblem(allencahn!, u0, (0.0, 0.5))

nsa() = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoFiniteDiff()))
inner_lincache(integ) = NonlinearSolveBase.get_linear_cache(integ.cache.nlsolver.cache.cache)

function reltol_after(alg, tol; nsteps = 5, kwargs...)
    integ = init(prob, alg; reltol = tol, abstol = tol, kwargs...)
    for _ in 1:nsteps
        step!(integ)
    end
    return inner_lincache(integ).reltol
end

@testset "the integrator tolerance reaches NonlinearSolveAlg's Krylov descent" begin
    # `build_nlsolver` zeroes the inner nonlinear tolerances and restores a fixed
    # `_inner_lintol` for the linear solve alone, so without the per-step update every
    # Krylov descent runs to ~1e-13 however loose the integrator is.
    for tol in (1.0e-3, 1.0e-6, 1.0e-9)
        @test reltol_after(FBDF(linsolve = KrylovJL_GMRES(), nlsolve = nsa()), tol) == tol
    end
    @test reltol_after(KenCarp4(linsolve = KrylovJL_GMRES(), nlsolve = nsa()), 1.0e-7) == 1.0e-7

    # A fixed-step run has no `reltol` to inherit; `NLNewton`'s `compute_step!` uses `eps`.
    @test reltol_after(
        FBDF(linsolve = KrylovJL_GMRES(), nlsolve = nsa()), 1.0e-6;
        adaptive = false, dt = 0.01
    ) == eps(Float64)
end

@testset "direct inner solves keep the fixed tolerance" begin
    # A factorization has no tolerance to set, and `update_tolerances!` throws for one.
    for ls in (LUFactorization(), nothing)
        alg = ls === nothing ? FBDF(nlsolve = nsa()) : FBDF(linsolve = ls, nlsolve = nsa())
        @test reltol_after(alg, 1.0e-3) == _inner_lintol(Float64)
    end
end

@testset "loosening the Krylov tolerance costs no accuracy" begin
    ref = solve(prob, Rodas5P(); reltol = 1.0e-13, abstol = 1.0e-15)
    for tol in (1.0e-6, 1.0e-9)
        sol = solve(
            prob, FBDF(linsolve = KrylovJL_GMRES(), nlsolve = nsa());
            reltol = tol, abstol = tol
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.stats.nnonlinconvfail == 0
        @test norm(sol.u[end] .- ref.u[end]) / norm(ref.u[end]) < 100tol
    end
end
