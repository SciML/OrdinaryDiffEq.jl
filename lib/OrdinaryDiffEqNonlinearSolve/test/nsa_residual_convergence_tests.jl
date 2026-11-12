using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using NonlinearSolve: DFSane, NewtonRaphson, Broyden, Klement, TrustRegion
using LineSearch: BackTracking
using ADTypes, LinearAlgebra, SciMLBase
using Test

# The outer convergence test used to infer "close to the root" from the displacement the inner
# solver produced, and never from the residual. That inference is only sound for a step that
# solves the Newton system, where the displacement is `W⁻¹F`. An inner solver whose line search
# collapses leaves the iterate unmoved — either stopping the cache with a failure retcode that
# `compute_step!` only looked for on *entry*, or returning a step orders of magnitude below
# roundoff with the residual bit-for-bit unchanged. Either way `ndz` is tiny for a reason that
# has nothing to do with the distance to the root, `iter == 1 && ndz < 1e-5` fires, and the
# stage is accepted unconverged.
#
# `DFSane` and a back-tracking `Broyden` genuinely cannot solve these stages. The property
# under test is only that the failure reaches the step controller instead of being dressed up
# as a converged solve: on unfixed code ROBER + FBDF + DFSane runs 21 accepted steps and
# returns the initial condition with `ReturnCode.Success`.

function rober!(du, u, p, t)
    y1, y2, y3 = u
    k1, k2, k3 = p
    du[1] = -k1 * y1 + k3 * y2 * y3
    du[2] = k1 * y1 - k2 * y2^2 - k3 * y2 * y3
    du[3] = k2 * y2^2
    return nothing
end
const rober_prob = ODEProblem(
    rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4)
)

function hires!(du, u, p, t)
    y1, y2, y3, y4, y5, y6, y7, y8 = u
    du[1] = -1.71y1 + 0.43y2 + 8.32y3 + 0.0007
    du[2] = 1.71y1 - 8.75y2
    du[3] = -10.03y3 + 0.43y4 + 0.035y5
    du[4] = 8.32y2 + 1.71y3 - 1.12y4
    du[5] = -1.745y5 + 0.43y6 + 0.43y7
    du[6] = -280.0y6 * y8 + 0.69y4 + 1.71y5 - 0.43y6 + 0.69y7
    du[7] = 280.0y6 * y8 - 1.81y7
    du[8] = -280.0y6 * y8 + 1.81y7
    return nothing
end
const hires_prob = ODEProblem(
    hires!, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0057], (0.0, 321.8122)
)

reference(prob) = solve(prob, Rodas5P(); reltol = 1.0e-13, abstol = 1.0e-15).u[end]
relerr(sol, ref) = norm(sol.u[end] .- ref) / norm(ref)

@testset "an unconverged stage is never reported as a successful solve" begin
    for (prob, ODEAlg, inner) in (
            (rober_prob, FBDF, () -> DFSane()),
            (rober_prob, KenCarp4, () -> DFSane()),
            (hires_prob, FBDF, () -> DFSane()),
            (hires_prob, FBDF, () -> Broyden(; linesearch = BackTracking())),
            (hires_prob, KenCarp4, () -> Broyden(; linesearch = BackTracking())),
        )
        ref = reference(prob)
        sol = solve(
            prob, ODEAlg(nlsolve = NonlinearSolveAlg(inner()));
            reltol = 1.0e-6, abstol = 1.0e-9, maxiters = 10^5
        )
        @test !SciMLBase.successful_retcode(sol) || relerr(sol, ref) < 1.0e-3
    end
end

# Quasi-Newton inners that cannot solve the stage must not report Success at `u0`;
# Newton-type inners solve it to `exp(-1)`.
@testset "quasi-Newton residual veto does not claim Success at u0" begin
    exact = exp(-1)
    prob = ODEProblem((du, u, p, t) -> (du .= -u), [1.0], (0.0, 1.0))
    for inner in (Broyden(), Klement())
        sol = solve(prob, TRBDF2(nlsolve = NonlinearSolveAlg(inner)))
        @test !SciMLBase.successful_retcode(sol)
        @test sol.u[end][1] ≈ 1.0 atol = 1.0e-6
    end
    for inner in (NewtonRaphson(), TrustRegion())
        sol = solve(prob, TRBDF2(nlsolve = NonlinearSolveAlg(inner)))
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end][1] ≈ exact rtol = 2.0e-2
    end
end

# `max_iter` bounds the inner steps of a stage and the iterations reported for it.
@testset "NonlinearSolveAlg.max_iter bounds inner steps and reported iterations" begin
    prob = ODEProblem((du, u, p, t) -> (du .= -u .^ 3), [1.0], (0.0, 1.0))
    for inner in (Broyden(), NewtonRaphson(), TrustRegion())
        integ = init(
            prob, ImplicitEuler(nlsolve = NonlinearSolveAlg(inner; max_iter = 1));
            dt = 1.0, adaptive = false, abstol = 1.0e-12, reltol = 1.0e-12
        )
        step!(integ)
        nc = integ.cache.nlsolver.cache.cache
        @test nc.nsteps <= 1
        @test integ.stats.nnonliniter <= 1
    end
end

# A back-tracking Broyden step moves the iterate by almost nothing on this stage, so
# the displacement test alone accepts every step at `u0` and the solve reports Success
# with `u = u0`. The residual veto makes the failure reach the step controller.
@testset "a collapsed quasi-Newton step is not a false Success at u0" begin
    prob = ODEProblem((du, u, p, t) -> (du .= -u), [1.0], (0.0, 1.0))
    sol = solve(
        prob,
        ImplicitEuler(nlsolve = NonlinearSolveAlg(Broyden(; linesearch = BackTracking())));
        maxiters = 100
    )
    @test !SciMLBase.successful_retcode(sol) ||
        isapprox(sol.u[end][1], exp(-1); rtol = 1.0e-2)
end

# Scaling the residual by a constant, through `precondition` or through a mass matrix
# with the right-hand side scaled to match, preserves every stage root. The veto must
# not read such a scaled residual against state tolerances.
@testset "root-preserving residual scaling keeps the stage root" begin
    prob = ODEProblem((u, p, t) -> -u .^ 2, [1.0], (0.0, 1.0))
    ref = 1.0
    for _ in 1:8
        ref = 2ref / (1 + sqrt(1 + 4 * (1 / 8) * ref))
    end
    for s in (1.0e-12, 1.0, 1.0e12, 1.0e24)
        G = let s = s
            (fu, u, p) -> s .* fu
        end
        sol = solve(
            prob,
            ImplicitEuler(nlsolve = NonlinearSolveAlg(TrustRegion(); precondition = G));
            dt = 1 / 8, adaptive = false, abstol = 1.0e-8, reltol = 1.0e-8
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end][1] ≈ ref atol = 1.0e-6
    end

    alg = TRBDF2(nlsolve = NonlinearSolveAlg(TrustRegion()))
    plain = ODEProblem((du, u, p, t) -> (du .= -u .^ 3), [1.0], (0.0, 1.0))
    scaled = ODEProblem(
        ODEFunction(
            (du, u, p, t) -> (du .= -1.0e18 .* u .^ 3);
            mass_matrix = Diagonal([1.0e18])
        ),
        [1.0], (0.0, 1.0)
    )
    sol_plain = solve(plain, alg; abstol = 1.0e-8, reltol = 1.0e-8)
    sol_scaled = solve(scaled, alg; abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol_scaled)
    @test sol_scaled.u[end][1] ≈ 1 / sqrt(3) atol = 1.0e-5
    @test sol_scaled.stats.nf <= 2 * sol_plain.stats.nf
end

# Each component of the residual is weighed against its own tolerance: a component
# with a loose tolerance does not lend it to a component with a tight one.
@testset "the residual veto weights components by their own tolerances" begin
    prob = ODEProblem((du, u, p, t) -> (du .= -u), [1.0, 1.0], (0.0, 1.0))
    integ = init(
        prob, ImplicitEuler(nlsolve = NonlinearSolveAlg(TrustRegion()));
        dt = 0.1, abstol = [1.0e-12, 1.0], reltol = 1.0e-10
    )
    OrdinaryDiffEqNonlinearSolve.initialize!(integ.cache.nlsolver, integ)
    nc = integ.cache.nlsolver.cache.cache
    OrdinaryDiffEqNonlinearSolve.get_u(nc) .= [0.0, 1.0]
    # `γΔt = 0.1` for this stage, so an `fu` entry of `1e-3` is a `1e-4` residual in
    # `z` units: far above the first component's tolerance, far below the second's.
    veto = OrdinaryDiffEqNonlinearSolve.unconverged_nlsolvealg_stage
    OrdinaryDiffEqNonlinearSolve.get_fu(nc) .= [1.0e-3, 0.0]
    @test veto(integ.cache.nlsolver, integ)
    OrdinaryDiffEqNonlinearSolve.get_fu(nc) .= [0.0, 1.0e-3]
    @test !veto(integ.cache.nlsolver, integ)
    OrdinaryDiffEqNonlinearSolve.get_fu(nc) .= [0.0, 0.0]
    @test !veto(integ.cache.nlsolver, integ)
end

@testset "a solvable stage still converges" begin
    for prob in (rober_prob, hires_prob), ODEAlg in (FBDF, KenCarp4)
        ref = reference(prob)
        sol = solve(
            prob,
            ODEAlg(nlsolve = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff())));
            reltol = 1.0e-6, abstol = 1.0e-9
        )
        @test SciMLBase.successful_retcode(sol)
        @test relerr(sol, ref) < 1.0e-4
    end
end
