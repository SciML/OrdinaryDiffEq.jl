using OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF
using NonlinearSolve, SciMLBase, Test
using ADTypes: AutoFiniteDiff, AutoForwardDiff

const NSA = OrdinaryDiffEqNonlinearSolve.NonlinearSolveAlg

# `u' = -k u` keeps the state inside `[exp(-1), 1]` over `tspan`, while the stage unknown
# `z` stays within one step of zero. The two ranges do not overlap, which is what lets the
# tests below tell which of the two a corrector is being handed.
scalar_prob(k = 1.0) = ODEProblem((u, p, t) -> -p.k * u, 1.0, (0.0, 1.0), (; k = k))
function vector_prob()
    return ODEProblem(
        (du, u, p, t) -> (du .= -p.k .* u), [1.0, 2.0], (0.0, 1.0), (; k = 1.0)
    )
end

@testset "corrector acts on the stage state, out-of-place" begin
    seen_u = Float64[]
    seen_uprev = Float64[]
    seen_p = Any[]
    function H(u, uprev, p, cache)
        push!(seen_u, u)
        push!(seen_uprev, uprev)
        push!(seen_p, p)
        return u
    end

    prob = scalar_prob()
    sol = solve(
        prob, ImplicitEuler(nlsolve = NSA(; postcondition = H));
        dt = 0.1, adaptive = false
    )

    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ exp(-1.0) rtol = 0.1
    @test !isempty(seen_u)
    # The state over this solve. `z` would be within 0.11 of zero throughout.
    @test all(x -> exp(-1.0) - 0.05 <= x <= 1.05, seen_u)
    # The ODE's parameters, not the stage problem's internal parameter tuple.
    @test all(==(prob.p), seen_p)

    # An identity corrector must reproduce the untouched solve exactly: the conjugation
    # back to `z` is written so that an unmodified state cancels bit-for-bit.
    sol_id = solve(
        prob, ImplicitEuler(nlsolve = NSA(; postcondition = (u, uprev, p, c) -> u));
        dt = 0.1, adaptive = false
    )
    sol_none = solve(prob, ImplicitEuler(nlsolve = NSA()); dt = 0.1, adaptive = false)
    @test sol_id.u == sol_none.u
end

@testset "corrector acts on the stage state, in-place" begin
    seen_u = Vector{Float64}[]
    seen_p = Any[]
    function H!(u, uprev, p, cache)
        push!(seen_u, copy(u))
        push!(seen_p, p)
        return nothing
    end

    prob = vector_prob()
    sol = solve(
        prob, ImplicitEuler(nlsolve = NSA(; postcondition = H!));
        dt = 0.1, adaptive = false
    )

    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ [exp(-1.0), 2exp(-1.0)] rtol = 0.1
    @test !isempty(seen_u)
    @test all(v -> exp(-1.0) - 0.05 <= v[1] <= 1.05, seen_u)
    @test all(v -> 2exp(-1.0) - 0.1 <= v[2] <= 2.1, seen_u)
    @test all(==(prob.p), seen_p)

    sol_id = solve(
        prob, ImplicitEuler(nlsolve = NSA(; postcondition = (u, uprev, p, c) -> nothing));
        dt = 0.1, adaptive = false
    )
    sol_none = solve(prob, ImplicitEuler(nlsolve = NSA()); dt = 0.1, adaptive = false)
    @test sol_id.u == sol_none.u
end

@testset "corrector acts on the stage state, multistep form" begin
    # BDF methods flip `nlsolver.method` to `COEFFICIENT_MULTISTEP` mid-solve, where the
    # stage unknown already *is* the state; the corrector must see the state either way.
    seen_u = Vector{Float64}[]
    function H!(u, uprev, p, cache)
        push!(seen_u, copy(u))
        return nothing
    end

    prob = vector_prob()
    sol = solve(
        prob, FBDF(nlsolve = NSA(; postcondition = H!)); dt = 0.05, adaptive = false
    )

    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ [exp(-1.0), 2exp(-1.0)] rtol = 0.05
    @test !isempty(seen_u)
    @test all(v -> exp(-1.0) - 0.05 <= v[1] <= 1.05, seen_u)
end

@testset "each stage of a multi-stage DIRK is seen at its own time point" begin
    # `u' = 1` from `u0 = 0` has exact state `t`, so the value the corrector is handed
    # identifies the time point of the stage it belongs to. TRBDF2's two implicit stages sit
    # at `t + γΔt` (γ = 2 - √2) and `t + Δt`.
    rec = Tuple{Bool, Float64}[]
    H(u, uprev, p, cache) = (push!(rec, (cache === nothing, u)); u)
    solve(
        ODEProblem((u, p, t) -> 1.0, 0.0, (0.0, 1.0)),
        TRBDF2(nlsolve = NSA(; postcondition = H)); dt = 0.5, adaptive = false
    )

    # The last correction of a stage's run sees that stage's accepted state; a run starts at
    # each predictor correction.
    finals = Float64[]
    for (i, (_, u)) in pairs(rec)
        (i == lastindex(rec) || rec[i + 1][1]) && push!(finals, u)
    end
    γ = 2 - sqrt(2)
    @test finals ≈ [0.5γ, 0.5, 0.5 + 0.5γ, 1.0]
end

@testset "predictor is corrected before the first residual" begin
    # The once-per-stage correction runs with no cache, and with `u == uprev`, which is the
    # stage analogue of the `H(u0, u0, p)` correction NonlinearSolve applies at `init`.
    firsts = Tuple{Float64, Float64}[]
    function H(u, uprev, p, cache)
        cache === nothing && push!(firsts, (u, uprev))
        return u
    end

    solve(
        scalar_prob(), ImplicitEuler(nlsolve = NSA(; postcondition = H));
        dt = 0.1, adaptive = false
    )
    @test length(firsts) == 10          # one per step
    @test all(t -> t[1] == t[2], firsts)
    @test all(t -> exp(-1.0) - 0.05 <= t[1] <= 1.05, firsts)
end

@testset "previous iterate of a stage's first correction is its predictor" begin
    # `reinit!` does not refresh the inner cache's previous-iterate buffer, so the raw
    # `u_prev` on the first commit of a stage is the last `z` of the *previous* stage read
    # through this stage's `tmp` — not a state at all (0.6 rather than 0.8 below).
    rec = Tuple{Bool, Float64, Float64}[]
    function H(u, uprev, p, cache)
        push!(rec, (cache === nothing, u, uprev))
        return u
    end

    solve(
        scalar_prob(), ImplicitEuler(nlsolve = NSA(; postcondition = H));
        dt = 0.25, adaptive = false
    )

    # Each stage's first in-loop correction must see that stage's corrected predictor.
    checked = 0
    predictor = nothing
    for (is_predictor, u, uprev) in rec
        if is_predictor
            predictor = u
        elseif predictor !== nothing
            @test uprev == predictor
            checked += 1
            predictor = nothing
        end
    end
    @test checked == 4      # one stage per fixed step
end

@testset "limiting corrector converges to the unlimited solution" begin
    # SPICE-style limiting: clip each accepted move relative to the previous iterate of the
    # same stage. The clip is inactive at a root, so the converged answer is unchanged.
    prob = scalar_prob(10.0)
    ref = solve(prob, ImplicitEuler(nlsolve = NSA()); abstol = 1.0e-6, reltol = 1.0e-6)

    function limited(δ, counter)
        return function (u, uprev, p, cache)
            du = u - uprev
            abs(du) > δ && (counter[] += 1)
            return uprev + clamp(du, -δ, δ)
        end
    end

    # Wide enough never to bind: the solve must then be the untouched one, bit for bit.
    loose = Ref(0)
    sol_loose = solve(
        prob, ImplicitEuler(nlsolve = NSA(; postcondition = limited(1.0e-2, loose)));
        abstol = 1.0e-6, reltol = 1.0e-6
    )
    @test loose[] == 0
    @test sol_loose.u == ref.u

    # Tight enough to clip repeatedly, and still converging to the same solution.
    tight = Ref(0)
    sol_tight = solve(
        prob, ImplicitEuler(nlsolve = NSA(; postcondition = limited(2.0e-3, tight)));
        abstol = 1.0e-6, reltol = 1.0e-6
    )
    @test SciMLBase.successful_retcode(sol_tight)
    @test tight[] > 0
    @test sol_tight.u[end] ≈ ref.u[end] rtol = 1.0e-8
end

@testset "precondition composes into the stage residual" begin
    # Row scaling is root-preserving and leaves the Newton iterates invariant, so it is the
    # transform whose effect on the solution is knowable in advance.
    prob = vector_prob()
    Gscale!(fu, u, p) = (fu[1] *= 2.0; fu[2] *= 3.0; nothing)
    sol = solve(
        prob, ImplicitEuler(nlsolve = NSA(; precondition = Gscale!));
        abstol = 1.0e-10, reltol = 1.0e-10
    )
    ref = solve(prob, ImplicitEuler(nlsolve = NSA()); abstol = 1.0e-10, reltol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ ref.u[end] rtol = 1.0e-10

    # Same coordinate contract as the corrector: stage state and the ODE's parameters.
    seen_u = Vector{Float64}[]
    seen_p = Any[]
    function Grec!(fu, u, p)
        push!(seen_u, copy(u))
        push!(seen_p, p)
        fu .*= 2.0
        return nothing
    end
    solr = solve(
        prob, ImplicitEuler(nlsolve = NSA(; precondition = Grec!)); dt = 0.1, adaptive = false
    )
    @test SciMLBase.successful_retcode(solr)
    @test !isempty(seen_u)
    @test all(v -> exp(-1.0) - 0.05 <= v[1] <= 1.05, seen_u)
    @test all(==(prob.p), seen_p)

    # Out-of-place form.
    sprob = scalar_prob()
    sols = solve(
        sprob, ImplicitEuler(nlsolve = NSA(; precondition = (fu, u, p) -> 2fu));
        abstol = 1.0e-10, reltol = 1.0e-10
    )
    sref = solve(sprob, ImplicitEuler(nlsolve = NSA()); abstol = 1.0e-10, reltol = 1.0e-10)
    @test sols.u[end] ≈ sref.u[end] rtol = 1.0e-10
end

@testset "precondition survives a no-init inner solver" begin
    # A SimpleNonlinearSolve inner algorithm lands in `NonlinearSolveNoInitCache`, which is
    # built twice — once with the integrator's zeroed tolerances and again without them —
    # and is driven by a complete `solve!` per outer iteration. The composed residual has to
    # survive both, and be composed exactly once.
    prob = vector_prob()
    inner = SimpleNewtonRaphson(autodiff = AutoFiniteDiff())
    seen_u = Vector{Float64}[]
    function G!(fu, u, p)
        push!(seen_u, copy(u))
        fu .*= 2.0
        return nothing
    end
    sol = solve(
        prob, ImplicitEuler(nlsolve = NSA(inner; precondition = G!));
        dt = 0.1, adaptive = false
    )
    # Row scaling is root-preserving, so the answer must match the unconditioned solve.
    ref = solve(prob, ImplicitEuler(nlsolve = NSA(inner)); dt = 0.1, adaptive = false)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ ref.u[end] rtol = 1.0e-8
    @test !isempty(seen_u)
    @test all(v -> exp(-1.0) - 0.05 <= v[1] <= 1.05, seen_u)
end

@testset "in-place precondition rejects dual-based inner autodiff" begin
    # Dropping `W` reuse leaves the inner solver differentiating the in-place stage
    # residual, which writes through preallocated `Float64` buffers. Better an explanation
    # up front than a `MethodError` from inside ForwardDiff.
    G!(fu, u, p) = (fu .*= 2.0; nothing)
    for inner in (
            NewtonRaphson(), NewtonRaphson(autodiff = AutoForwardDiff()),
            SimpleNewtonRaphson(), TrustRegion(),
        )
        @test_throws ArgumentError init(
            vector_prob(), ImplicitEuler(nlsolve = NSA(inner; precondition = G!)); dt = 0.1
        )
    end
    # Only `precondition` disables `W` reuse, so a corrector is unaffected by the backend.
    for inner in (NewtonRaphson(), NewtonRaphson(autodiff = AutoForwardDiff()), TrustRegion())
        @test init(
            vector_prob(),
            ImplicitEuler(nlsolve = NSA(inner; postcondition = (u, up, p, c) -> nothing));
            dt = 0.1
        ) isa SciMLBase.DEIntegrator
    end
    # Out-of-place residuals allocate, so every backend works there.
    sol = solve(
        scalar_prob(),
        ImplicitEuler(nlsolve = NSA(NewtonRaphson(); precondition = (fu, u, p) -> 2fu));
        dt = 0.1, adaptive = false
    )
    @test SciMLBase.successful_retcode(sol)
end

@testset "precondition disables W reuse" begin
    prob = vector_prob()
    plain = init(prob, ImplicitEuler(nlsolve = NSA()); dt = 0.1, adaptive = false)
    @test plain.cache.nlsolver.cache.W !== nothing

    G!(fu, u, p) = (fu .*= 2.0; nothing)
    cond = init(
        prob, ImplicitEuler(nlsolve = NSA(; precondition = G!)); dt = 0.1, adaptive = false
    )
    # `W` is the Jacobian of the raw stage residual, not of the composed one.
    @test cond.cache.nlsolver.cache.W === nothing

    post = init(
        prob, ImplicitEuler(nlsolve = NSA(; postcondition = (u, up, p, c) -> nothing));
        dt = 0.1, adaptive = false
    )
    # A corrector does not change the residual, so `W` reuse stays valid.
    @test post.cache.nlsolver.cache.W !== nothing
end

@testset "corrector survives a resize" begin
    seen_lengths = Int[]
    function H!(u, uprev, p, cache)
        push!(seen_lengths, length(u))
        return nothing
    end
    function grow!(integ)
        resize!(integ, length(integ.u) + 1)
        integ.u[end] = 1.0
        return nothing
    end
    cb = DiscreteCallback((u, t, integ) -> t >= 0.5 && length(integ.u) == 2, grow!)

    sol = solve(
        vector_prob(), ImplicitEuler(nlsolve = NSA(; postcondition = H!));
        dt = 0.1, adaptive = false, callback = cb
    )
    @test SciMLBase.successful_retcode(sol)
    @test length(sol.u[end]) == 3
    @test sort(unique(seen_lengths)) == [2, 3]
end

@testset "unsupported nonlinear solvers reject rather than ignore" begin
    H(u, uprev, p, cache) = u
    G(fu, u, p) = fu
    for T in (
            OrdinaryDiffEqNonlinearSolve.NLNewton,
            OrdinaryDiffEqNonlinearSolve.NLFunctional,
            OrdinaryDiffEqNonlinearSolve.NLAnderson,
            OrdinaryDiffEqNonlinearSolve.HomotopyNonlinearSolveAlg,
        )
        @test_throws ArgumentError T(postcondition = H)
        @test_throws ArgumentError T(precondition = G)
        @test T() isa OrdinaryDiffEqNonlinearSolve.AbstractNLSolverAlgorithm
    end

    # A NonlinearSolve.jl algorithm that does not apply iterate corrections is rejected by
    # NonlinearSolve's own `unsupported_postcondition` report when the stage cache is built.
    @test_throws Exception init(
        vector_prob(),
        ImplicitEuler(nlsolve = NSA(SimpleNewtonRaphson(); postcondition = H));
        dt = 0.1
    )
end

@testset "DAE stage systems reject conditioning" begin
    H(u, uprev, p, cache) = u
    daef(out, du, u, p, t) = (out[1] = du[1] + u[1]; nothing)
    prob = DAEProblem(daef, [-1.0], [1.0], (0.0, 1.0); differential_vars = [true])
    @test_throws ArgumentError init(
        prob, DImplicitEuler(nlsolve = NSA(; postcondition = H)); dt = 0.1
    )
end
