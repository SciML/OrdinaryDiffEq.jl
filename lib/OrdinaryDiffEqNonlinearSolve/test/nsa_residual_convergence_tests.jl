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
        @test !SciMLBase.successful_retcode(sol) ||
            isapprox(sol.u[end][1], exact; rtol = 2.0e-2)
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
    # A residual-based trust radius shrinks with `s`, making an unsolved stage look
    # stationary at small scales. Fix the radius in state units to test the veto.
    inner = TrustRegion(max_trust_radius = 1.0, initial_trust_radius = 1.0)
    for s in (1.0e-12, 1.0, 1.0e12, 1.0e24)
        G = let s = s
            (fu, u, p) -> s .* fu
        end
        sol = solve(
            prob,
            ImplicitEuler(nlsolve = NonlinearSolveAlg(inner; precondition = G));
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

# A non-finite residual must not read as convergence through `Inf <= eps * Inf` or `NaN > 1`.
@testset "a non-finite residual is not a solved stage" begin
    prob = ODEProblem((du, u, p, t) -> (du .= -u), [1.0, 1.0], (0.0, 1.0))
    integ = init(
        prob, ImplicitEuler(nlsolve = NonlinearSolveAlg(TrustRegion()));
        dt = 0.1, abstol = 1.0e-6, reltol = 1.0e-3
    )
    OrdinaryDiffEqNonlinearSolve.initialize!(integ.cache.nlsolver, integ)
    nc = integ.cache.nlsolver.cache.cache
    OrdinaryDiffEqNonlinearSolve.get_u(nc) .= [0.0, 1.0]
    veto = OrdinaryDiffEqNonlinearSolve.unconverged_nlsolvealg_stage
    for bad in (Inf, -Inf, NaN)
        OrdinaryDiffEqNonlinearSolve.get_fu(nc) .= [bad, 0.0]
        @test veto(integ.cache.nlsolver, integ)
    end

    # One implicit-Euler step of `u' = -u`, `u0 = 1`, `dt = 1` has stage root `1 / 2`.
    for bad in (Inf, NaN)
        f! = let bad = bad
            (du, u, p, t) -> (du[1] = u[1] > 1 ? bad : -u[1]; nothing)
        end
        prob = ODEProblem(f!, [1.0], (0.0, 1.0))
        sol = solve(
            prob,
            ImplicitEuler(nlsolve = NonlinearSolveAlg(Broyden(alpha = 1.0e12)));
            dt = 1.0, adaptive = false, abstol = 1.0e-6, reltol = 1.0e-3
        )
        @test !SciMLBase.successful_retcode(sol) ||
            isapprox(sol.u[end][1], 0.5; atol = 1.0e-3)
    end
end

@testset "finite stage terms do not overflow the roundoff floor" begin
    residual_or_zero = OrdinaryDiffEqNonlinearSolve._stage_resid_or_zero
    for T in (Float32, Float64)
        z = floatmax(T) * T(0.75)
        for r in (z, -z)
            @test residual_or_zero(r, r) == r
            @test residual_or_zero(r, -r) == r
            @test iszero(residual_or_zero(eps(T) * r, r))
        end
    end

    for inplace in (false, true)
        f = inplace ? (du, u, p, t) -> (du .= -u) : (u, p, t) -> -u
        prob = ODEProblem{inplace}(f, [1.0], (0.0, 1.0))
        integ = init(
            prob, ImplicitEuler(nlsolve = NonlinearSolveAlg(TrustRegion()));
            dt = 1.0, abstol = 1.0e-6, reltol = 1.0e-3
        )
        OrdinaryDiffEqNonlinearSolve.initialize!(integ.cache.nlsolver, integ)
        nc = integ.cache.nlsolver.cache.cache
        OrdinaryDiffEqNonlinearSolve.get_u(nc) .= 1.0e308
        OrdinaryDiffEqNonlinearSolve.get_fu(nc) .= 1.0e308
        @test OrdinaryDiffEqNonlinearSolve.unconverged_nlsolvealg_stage(
            integ.cache.nlsolver, integ
        )
    end
end

# A one-ulp change of the stage state moves `γΔt · f` by about `γΔt · λ · eps`, far above
# the tolerances, so a stiff stage solved to roundoff still has a large residual.
@testset "a stiff stage at its roundoff floor is converged" begin
    λ, h = 1.0e14, 1 / 8
    exact = 1 - h / (1 + h * λ)
    rhs_iip = (du, u, p, t) -> (du[1] = λ * (1 - u[1]) - 1; nothing)
    rhs_oop = (u, p, t) -> [λ * (1 - u[1]) - 1]
    inners = (TrustRegion(max_trust_radius = 1.0, initial_trust_radius = 1.0), Broyden())
    for rhs in (rhs_iip, rhs_oop), inner in inners
        sol = solve(
            ODEProblem(rhs, [1.0], (0.0, h)), ImplicitEuler(nlsolve = NonlinearSolveAlg(inner));
            dt = h, adaptive = false, abstol = 1.0e-8, reltol = 1.0e-8
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end][1] ≈ exact atol = 1.0e-8
    end
    for rhs in (rhs_iip, rhs_oop), ODEAlg in (TRBDF2, FBDF)
        sol = solve(
            ODEProblem(rhs, [1.0], (0.0, h)), ODEAlg(nlsolve = NonlinearSolveAlg(TrustRegion()));
            abstol = 1.0e-8, reltol = 1.0e-8
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end][1] ≈ 1 - 1 / λ atol = 1.0e-8
        @test sol.stats.nf <= 1000
    end

    # A smooth state shift nearly cancels in a discrete Laplacian, whose roundoff floor
    # is set by its largest entries.
    N = 2000
    dx = 1 / (N + 1)
    lap! = (du, u, p, t) -> begin
        for i in 1:N
            l = i > 1 ? u[i - 1] : 0.0
            r = i < N ? u[i + 1] : 0.0
            du[i] = (l - 2u[i] + r) / dx^2
        end
        nothing
    end
    u0 = [sin(pi * i * dx) + 1.0 for i in 1:N]
    A = SymTridiagonal(fill(-2.0, N), ones(N - 1)) / dx^2
    sol = solve(
        ODEProblem(lap!, u0, (0.0, 1.0)),
        ImplicitEuler(nlsolve = NonlinearSolveAlg(TrustRegion()));
        dt = 1.0, adaptive = false, abstol = 1.0e-10, reltol = 1.0e-10
    )
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ (I - A) \ u0 atol = 1.0e-8
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
