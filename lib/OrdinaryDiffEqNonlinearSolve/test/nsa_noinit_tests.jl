# NonlinearSolveAlg with inner algorithms that have no `init` (every SimpleNonlinearSolve
# algorithm): these land in `NonlinearSolveBase.NonlinearSolveNoInitCache`, a fallback with
# no iteration state, so `step!`, `get_fu`, `.stats` and `not_terminated` are all
# unavailable and the cache can only be driven by complete `solve!` calls. A
# SimpleNewtonRaphson solve ends on the integrator's own κ/η test; the others end on their
# own default criterion.
using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using NonlinearSolve: NewtonRaphson
using SimpleNonlinearSolve
using NonlinearSolveBase
using SciMLBase, LinearAlgebra
using Test

prob_scalar = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))

function vdp!(du, u, p, t)
    du[1] = u[2]
    du[2] = p * ((1 - u[1]^2) * u[2]) - u[1]
    return nothing
end
prob_vdp = ODEProblem(vdp!, [2.0, 0.0], (0.0, 1.0), 10.0)

function rober_mm!(du, u, p, t)
    y₁, y₂, y₃ = u
    du[1] = -0.04y₁ + 1.0e4 * y₂ * y₃
    du[2] = 0.04y₁ - 1.0e4 * y₂ * y₃ - 3.0e7 * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    return nothing
end
prob_mm = ODEProblem(
    ODEFunction(rober_mm!, mass_matrix = Diagonal([1.0, 1.0, 0.0])),
    [1.0, 0.0, 0.0], (0.0, 0.1)
)

noinit_cache(integrator) = integrator.cache.nlsolver.cache

@testset "cache-type identification" begin
    integ = init(prob_vdp, TRBDF2(nlsolve = NonlinearSolveAlg(SimpleNewtonRaphson())))
    nlc = noinit_cache(integ)
    @test nlc.cache isa NonlinearSolveBase.NonlinearSolveNoInitCache
    # W reuse stays enabled for no-init caches in place: the reused W is handed to the
    # inner solver as `prob.f.jac`, which SimpleNonlinearSolve honours, and the residual
    # writes through preallocated Float64 buffers that AD cannot go through.
    @test nlc.W !== nothing
    integ_oop = init(prob_scalar, TRBDF2(nlsolve = NonlinearSolveAlg(SimpleNewtonRaphson())))
    nlc_oop = noinit_cache(integ_oop)
    @test nlc_oop.cache isa NonlinearSolveBase.NonlinearSolveNoInitCache
    @test nlc_oop.W === nothing
    integ_tr = init(prob_vdp, TRBDF2(nlsolve = NonlinearSolveAlg(SimpleTrustRegion())))
    @test noinit_cache(integ_tr).cache isa NonlinearSolveBase.NonlinearSolveNoInitCache
    integ_nr = init(prob_vdp, TRBDF2(nlsolve = NonlinearSolveAlg(NewtonRaphson())))
    @test !(noinit_cache(integ_nr).cache isa NonlinearSolveBase.NonlinearSolveNoInitCache)
end

inner_algs = [
    ("SimpleNewtonRaphson", SimpleNewtonRaphson()),
    ("SimpleTrustRegion", SimpleTrustRegion()),
    ("SimpleKlement", SimpleKlement()),
    ("SimpleBroyden", SimpleBroyden()),
]

@testset "functional matrix: $pname" for (pname, prob, reltol) in [
        ("scalar oop", prob_scalar, 1.0e-5),
        ("vdp iip", prob_vdp, 1.0e-6),
        ("rober-mm iip", prob_mm, 1.0e-7),
    ]
    ref = solve(prob, Rodas5P(); reltol = 1.0e-12, abstol = 1.0e-12)
    for method in (Trapezoid, TRBDF2, FBDF), (iname, ialg) in inner_algs
        # SimpleBroyden diverges on this one cell (a clear failure, not a wrong answer):
        # its secant update cannot hold the algebraic row of the singular mass matrix.
        pname == "rober-mm iip" && method === Trapezoid && iname == "SimpleBroyden" &&
            continue
        sol = solve(
            prob, method(nlsolve = NonlinearSolveAlg(ialg));
            reltol = 1.0e-8, abstol = 1.0e-10
        )
        @test SciMLBase.successful_retcode(sol.retcode)
        refu = ref(sol.t[end])
        @test norm(sol.u[end] .- refu) / max(norm(refu), 1.0e-16) < reltol
    end
end

@testset "tolerances" begin
    # The build paths zero the inner tolerances so the integrator owns convergence on the
    # `step!`-driven path — but a no-init cache is `solve!`-driven, so zeroed tolerances
    # would leave every complete inner solve returning `MaxIters` (except when the residual
    # lands on exactly 0.0, which is why easy problems cannot see the bug). The rebuilt
    # cache must carry neither the zeroed tolerances nor the stepped-path-only kwargs that
    # would be splatted into `solve`.
    for integ in (
            init(prob_vdp, TRBDF2(nlsolve = NonlinearSolveAlg(SimpleNewtonRaphson()))),
            init(prob_scalar, TRBDF2(nlsolve = NonlinearSolveAlg(SimpleNewtonRaphson()))),
        )
        nlcache = noinit_cache(integ).cache
        @test !iszero(NonlinearSolveBase.get_abstol(nlcache))
        @test !iszero(NonlinearSolveBase.get_reltol(nlcache))
        @test !haskey(nlcache.kwargs, :termination_condition)
        @test !haskey(nlcache.kwargs, :linsolve_kwargs)
    end
end

@testset "resize" begin
    condition(u, t, integrator) = t == 0.5
    function affect!(integrator)
        resize!(integrator, 2)
        integrator.u[2] = 1.0
        return nothing
    end
    cb = DiscreteCallback(condition, affect!)
    decay!(du, u, p, t) = (@. du = -u; nothing)
    prob_grow = ODEProblem(decay!, [1.0], (0.0, 1.0))
    for (iname, ialg) in [
            ("SimpleNewtonRaphson", SimpleNewtonRaphson()),
            ("SimpleTrustRegion", SimpleTrustRegion()),
        ]
        integ = init(
            prob_grow, TRBDF2(nlsolve = NonlinearSolveAlg(ialg));
            callback = cb, tstops = [0.5], reltol = 1.0e-8, abstol = 1.0e-10
        )
        sol = solve!(integ)
        @test SciMLBase.successful_retcode(sol.retcode)
        # The second component decays from 1.0 at t = 0.5; a stale 1x1 `WReuseJac` left
        # aliased across the resize surfaces as a `SingularException` here (and would be
        # a silently wrong Jacobian at other sizes).
        @test isapprox(sol.u[end][1], exp(-1); atol = 1.0e-5)
        @test isapprox(sol.u[end][2], exp(-0.5); atol = 1.0e-5)
        nlcache = noinit_cache(integ).cache
        @test nlcache isa NonlinearSolveBase.NonlinearSolveNoInitCache
        # The post-resize rebuild must preserve the no-init tolerances exactly as the
        # build path does; the only termination criterion is the per-solve one.
        @test !iszero(NonlinearSolveBase.get_abstol(nlcache))
        if ialg isa SimpleNewtonRaphson
            @test nlcache.kwargs[:termination_condition] isa
                OrdinaryDiffEqNonlinearSolve.StageConvergenceMode
        else
            @test !haskey(nlcache.kwargs, :termination_condition)
        end
    end
end

@testset "rejected trial step" begin
    integ = init(
        prob_vdp, TRBDF2(nlsolve = NonlinearSolveAlg(SimpleTrustRegion()));
        reltol = 1.0e-8, abstol = 1.0e-10
    )
    step!(integ)
    nls = integ.cache.nlsolver
    @test nls.cache.cache isa NonlinearSolveBase.NonlinearSolveNoInitCache
    # A no-init cache has no trial-step state to reject: a zero displacement after a
    # complete `solve!` is genuine convergence, and `not_terminated` has nothing to read.
    @test !OrdinaryDiffEqNonlinearSolve._uninformative_step(nls, 0.0)
    @test !OrdinaryDiffEqNonlinearSolve._uninformative_step(nls, 1.0)
    # The stepped path keeps its #3817 discrimination untouched.
    integ_nr = init(
        prob_vdp, TRBDF2(nlsolve = NonlinearSolveAlg(NewtonRaphson()));
        reltol = 1.0e-8, abstol = 1.0e-10
    )
    step!(integ_nr)
    nls_nr = integ_nr.cache.nlsolver
    @test OrdinaryDiffEqNonlinearSolve._uninformative_step(nls_nr, 0.0) ==
        (NonlinearSolveBase.not_terminated(nls_nr.cache.cache) || nls_nr.cache.stalled)
    @test !OrdinaryDiffEqNonlinearSolve._uninformative_step(nls_nr, 1.0)
    sol = solve(
        prob_vdp, TRBDF2(nlsolve = NonlinearSolveAlg(SimpleTrustRegion()));
        reltol = 1.0e-8, abstol = 1.0e-10
    )
    @test SciMLBase.successful_retcode(sol.retcode)
end

# ROBER in units scaled by `s` (`y = s⋅Y`): identical dynamics, so tolerances scaled with it
# must give the same steps and the same weighted stage accuracy at every `s`.
function rober_scaled(u, s)
    y₁, y₂, y₃ = u ./ s
    return s .* [
        -0.04y₁ + 1.0e4 * y₂ * y₃,
        0.04y₁ - 1.0e4 * y₂ * y₃ - 3.0e7 * y₂^2,
        3.0e7 * y₂^2,
    ]
end
function rober_scaled_jac(u, s)
    y₁, y₂, y₃ = u ./ s
    return [
        -0.04 1.0e4 * y₃ 1.0e4 * y₂
        0.04 (-1.0e4 * y₃ - 6.0e7 * y₂) -1.0e4 * y₂
        0.0 6.0e7 * y₂ 0.0
    ]
end
rober_scaled!(du, u, s, t) = (du .= rober_scaled(u, s); nothing)
rober_scaled_oop(u, s, t) = rober_scaled(u, s)

# Weighted distance of an accepted ImplicitEuler step from the exact root of
# `u = uprev + h⋅f(u)`, in the integrator's own norm.
function implicit_euler_stage_error(integ, s)
    uprev, u, h = collect(integ.uprev), collect(integ.u), integ.t - integ.tprev
    x = copy(u)
    for _ in 1:30
        x = x .- (I - h * rober_scaled_jac(x, s)) \ (x .- uprev .- h .* rober_scaled(x, s))
    end
    w = integ.opts.abstol .+ integ.opts.reltol .* max.(abs.(uprev), abs.(u))
    return sqrt(sum(abs2, (u .- x) ./ w) / length(u))
end

# The integrator's κ/η test accepts a stage once its estimated weighted error is below
# κ = 1/100 (NLNewton stays under 0.04 here), so an accepted stage a whole tolerance unit off
# the root was accepted by some other criterion.
@testset "inner termination follows the integrator's tolerances (s = $s)" for s in
    (1.0e-10, 1.0e6)
    tspan = (0.0, 1.0e3)
    reltol, abstol = 1.0e-4, 1.0e-8 * s
    ref_steps = length(
        solve(
            ODEProblem(rober_scaled!, [s, 0.0, 0.0], tspan, s), ImplicitEuler();
            reltol, abstol
        ).t
    )
    for iip in (true, false)
        prob = iip ? ODEProblem(rober_scaled!, [s, 0.0, 0.0], tspan, s) :
            ODEProblem(rober_scaled_oop, [s, 0.0, 0.0], tspan, s)
        integ = init(
            prob, ImplicitEuler(nlsolve = NonlinearSolveAlg(SimpleNewtonRaphson()));
            reltol, abstol, maxiters = 3 * ref_steps
        )
        worst = 0.0
        for _ in integ
            # The iterator also yields the state a `MaxIters` abort leaves behind, whose `u` is
            # a rejected trial value rather than a stage the nonlinear solver accepted.
            integ.accept_step || continue
            worst = max(worst, implicit_euler_stage_error(integ, s))
        end
        @testset "$(iip ? "iip" : "oop")" begin
            @test worst < 1
            @test SciMLBase.successful_retcode(integ.sol.retcode)
        end
    end
end

# Weighted RMS distance of one fixed ImplicitEuler step from the stage root.
function one_step_stage_error(f, f!, u0, root, iip, inner; abstol, reltol)
    prob = ODEProblem(iip ? f! : f, u0, (0.0, 1.0))
    integ = init(
        prob, ImplicitEuler(nlsolve = NonlinearSolveAlg(inner));
        dt = 1.0, adaptive = false, abstol, reltol
    )
    step!(integ)
    w = abstol .+ reltol .* max.(abs.(integ.uprev), abs.(integ.u))
    return integ.sol.retcode, sqrt(sum(abs2, (integ.u .- root) ./ w) / length(root))
end

# One large component must not set the convergence floor for a small one: the stage root of
# `u' = [-u₁, -u₂²]`, `u₀ = [1e16, 1]`, `h = 1` is `[5e15, (√5 - 1)/2]`, and the Jacobian is
# well conditioned.
@testset "mixed-scale stage ($(nameof(typeof(inner))) $(iip ? "iip" : "oop"))" for inner in
        (SimpleNewtonRaphson(), SimpleTrustRegion()), iip in (true, false)
    f(u, p, t) = [-u[1], -u[2]^2]
    f!(du, u, p, t) = (du .= f(u, p, t); nothing)
    retcode, err = one_step_stage_error(
        f, f!, [1.0e16, 1.0], [5.0e15, (sqrt(5.0) - 1) / 2], iip, inner;
        abstol = 1.0e-8, reltol = 1.0e-6
    )
    @test SciMLBase.successful_retcode(retcode)
    @test err < 1
end

# A residual at the roundoff floor says nothing when `W` is ill-conditioned: the stage of
# `u' = [(1 - δ)u₁ + 2δ, 1]` from `u₀ = 0` with `h = 1` is linear, `W = diag(δ, 1)`, and its root
# is exactly `[2, 1]`.
@testset "ill-conditioned stage (δ = $δ, $(iip ? "iip" : "oop"))" for δ in
        (1.0e-12, eps(1.0)), iip in (true, false)
    f(u, p, t) = [(1 - δ) * u[1] + 2δ, 1.0]
    f!(du, u, p, t) = (du .= f(u, p, t); nothing)
    root = [2δ / (1 - (1 - δ)), 1.0]
    retcode, err = one_step_stage_error(
        f, f!, [0.0, 0.0], root, iip, SimpleNewtonRaphson();
        abstol = 1.0e-10, reltol = 1.0e-8
    )
    @test SciMLBase.successful_retcode(retcode)
    @test err < 1
end

# A stage whose Newton correction rounds to zero is converged, not stuck: warm-started at the
# previous inner solve's root, the stage of `u' = -u²` from `u₀ = 1`, `h = 1`
# (root `(√5 - 1)/2`) must be accepted down to `tol = eps`.
@testset "zero first correction ($(iip ? "iip" : "oop"), tol = $tol)" for iip in (true, false),
        tol in (1.0e-12, eps(1.0))
    f(u, p, t) = -u .^ 2
    f!(du, u, p, t) = (du .= f(u, p, t); nothing)
    integ = init(
        ODEProblem(iip ? f! : f, [1.0], (0.0, 1.0)),
        ImplicitEuler(nlsolve = NonlinearSolveAlg(SimpleNewtonRaphson()));
        dt = 1.0, adaptive = false, abstol = tol, reltol = tol
    )
    step!(integ)
    @test SciMLBase.successful_retcode(integ.sol.retcode)
    @test integ.t == 1.0
    @test integ.u[1] ≈ (sqrt(5.0) - 1) / 2 rtol = 1.0e-12
end

# `sin(x) = 0` from floating-point π: the Newton correction is stationary at rounding level,
# so only a first-iteration certificate can end the solve.
mutable struct ZeroCorrectionOwner
    inner_ndzprev::Float64
end
@testset "stationary roundoff" begin
    owner = ZeroCorrectionOwner(-1.0)
    norm_w = OrdinaryDiffEqNonlinearSolve.StageDisplacementNorm(
        [3.0], 1.0e-12, 1.0e-12, (x, t) -> norm(x), 0.0
    )
    newton_step(fu, x) = fu ./ cos.(x)
    mode = OrdinaryDiffEqNonlinearSolve.StageConvergenceMode(norm_w, 0.01, owner, newton_step)
    sol = solve(
        NonlinearProblem((x, p) -> sin.(x), [Float64(pi)]), SimpleNewtonRaphson();
        termination_condition = mode
    )
    @test SciMLBase.successful_retcode(sol.retcode)
    @test abs(sol.u[1] - Float64(pi)) <= eps(Float64(pi))
end

# A correction that rounds away is not a converged one. The stage of
# `u' = (u - 10¹⁶) - (u - 2) - 10⁻⁸(u - 2)²` from `u₀ = 10¹⁶`, `h = 1`, has the exact root 2, but
# the predictor 4 sits where the in-place reused-`W` correction is below half an ulp of the
# increment `z ≈ -10¹⁶`: the stored iterate does not move. The stage must not be accepted at 4.
@testset "rounded-away correction ($(iip ? "iip" : "oop"))" for iip in (true, false)
    function f(u, p, t)
        x = u[1] - 2.0
        return [(u[1] - 1.0e16) - (x + 1.0e-8 * x * x)]
    end
    f!(du, u, p, t) = (du .= f(u, p, t); nothing)
    predictor = iip ? ((u, uprev, p, t, dt) -> (u[1] = 4.0; nothing)) :
        ((uprev, p, t, dt) -> [4.0])
    integ = init(
        ODEProblem(iip ? f! : f, [1.0e16], (0.0, 1.0)),
        ImplicitEuler(nlsolve = NonlinearSolveAlg(SimpleNewtonRaphson(); predictor));
        dt = 1.0, adaptive = false, abstol = 1.0e-12, reltol = 0.0
    )
    step!(integ)
    @test !SciMLBase.successful_retcode(integ.sol.retcode) ||
        abs(integ.u[1] - 2.0) <= 1.0e-12
end
