using OrdinaryDiffEq, OrdinaryDiffEqSDIRK, Test
using OrdinaryDiffEqBDF: FBDFCache
using OrdinaryDiffEqCore: default_subcaches, full_cache
using ADTypes: AutoForwardDiff, AutoSparse
using SparseConnectivityTracer: TracerLocalSparsityDetector
using SparseMatrixColorings: GreedyColoringAlgorithm

decay!(du, u, p, t) = (du .= -p .* u; nothing)
decay_prefix!(du, u, p, t) = (du .= -view(p, 1:length(u)) .* u; nothing)

const DELETEAT_ADDAT_ALGS = (
    Trapezoid(), ImplicitEuler(), KenCarp4(), Tsit5(), Vern7(), Rodas5P(),
)

function reduced_reference(integ, alg, p, tend; kwargs...)
    prob = ODEProblem(decay!, copy(integ.u), (integ.t, tend), p)
    return solve(prob, alg; kwargs...).u[end]
end

function step_history(integ)
    uprev3 = hasfield(typeof(integ.cache), :uprev3) ? integ.cache.uprev3 : nothing
    return filter(v -> v isa AbstractVector && v !== integ.u, Any[integ.uprev, integ.uprev2, uprev3])
end

@testset "deleteat!/addat! mid-solve with $(nameof(typeof(alg)))" for alg in
    DELETEAT_ADDAT_ALGS
    tend = 1.0
    kwargs = (; adaptive = false, dt = 1 / 64, save_everystep = false)

    @testset "deleteat!" begin
        p = collect(1.0:6.0)
        integ = init(ODEProblem(decay!, collect(1.0:6.0), (0.0, tend), p), alg; kwargs...)
        for _ in 1:5
            step!(integ)
        end
        kept = integ.u[[1, 4, 5, 6]]
        deleteat!(integ, [2, 3])
        deleteat!(p, [2, 3])
        @test integ.u == kept
        ref = reduced_reference(integ, alg, copy(p), tend; kwargs...)
        solve!(integ)
        @test integ.sol.retcode == ReturnCode.Success
        @test integ.t == tend
        @test integ.u ≈ ref rtol = 1.0e-12
    end

    @testset "addat!" begin
        p = collect(1.0:6.0)
        integ = init(ODEProblem(decay!, collect(1.0:6.0), (0.0, tend), p), alg; kwargs...)
        for _ in 1:5
            step!(integ)
        end
        old = copy(integ.u)
        addat!(integ, 7:8)
        append!(p, [7.0, 8.0])
        @test integ.u[1:6] == old
        integ.u[7:8] .= [10.0, 20.0]
        u_modified!(integ, true)
        ref = reduced_reference(integ, alg, copy(p), tend; kwargs...)
        solve!(integ)
        @test integ.sol.retcode == ReturnCode.Success
        @test integ.t == tend
        @test integ.u ≈ ref rtol = 1.0e-12
    end

    @testset "addat! with a value" begin
        p = collect(1.0:6.0)
        integ = init(ODEProblem(decay!, collect(1.0:6.0), (0.0, tend), p), alg; kwargs...)
        for _ in 1:5
            step!(integ)
        end
        old = copy(integ.u)
        addat!(integ, 7:8, [10.0, 20.0])
        append!(p, [7.0, 8.0])
        @test integ.u == [old; 10.0; 20.0]
        @test all(v -> v[7:8] == [10.0, 20.0], step_history(integ))
        ref = reduced_reference(integ, alg, copy(p), tend; kwargs...)
        solve!(integ)
        @test integ.sol.retcode == ReturnCode.Success
        @test integ.u ≈ ref rtol = 1.0e-12
    end

    @testset "addat! without a value starts the new components at zero" begin
        integ = init(ODEProblem(decay!, collect(1.0:6.0), (0.0, tend), collect(1.0:6.0)), alg; kwargs...)
        for _ in 1:5
            step!(integ)
        end
        # Shrinking keeps a Vector's contents, so without a fill the old values come back.
        deleteat!(integ, 5:6)
        addat!(integ, 5:6)
        @test iszero(integ.u[5:6])
        @test all(v -> iszero(v[5:6]), step_history(integ))
        solve!(integ)
        @test integ.sol.retcode == ReturnCode.Success
        @test iszero(integ.u[5:6])
    end
end

@testset "addat! only appends" begin
    integ = init(ODEProblem(decay!, [1.0, 2.0], (0.0, 1.0), [1.0, 2.0]), Tsit5())
    step!(integ)
    @test_throws ArgumentError addat!(integ, 1:1)
    @test_throws ArgumentError addat!(integ, [3, 5])
    @test_throws DimensionMismatch addat!(integ, 3:4, [1.0, 2.0, 3.0])
    @test_throws DimensionMismatch addat!(integ, 3:4, [1.0 2.0])
    @test_throws InexactError addat!(integ, 3:3, 1.0im)
    @test length(integ.u) == 2
    addat!(integ, Int[])
    @test length(integ.u) == 2
    addat!(integ, 3)
    @test length(integ.u) == 3
end

@testset "default solver switches to a stiff method after $name" for (name, grow!) in (
        "addat!" => integ -> addat!(integ, 7:8),
        "resize!" => integ -> (resize!(integ, 8); integ.u[7:8] .= 0.0),
    )
    prob = ODEProblem(decay_prefix!, collect(1.0:6.0), (0.0, 1.0), collect(1.0:8.0))
    cb = DiscreteCallback((u, t, integ) -> t >= 0.3 && length(u) == 6, grow!)
    sol = solve(prob; callback = cb, abstol = 1.0e-8, reltol = 1.0e-8)
    @test sol.retcode == ReturnCode.Success
    @test 4 in sol.alg_choice
    @test iszero(sol.u[end][7:8])
end

@testset "addat! writes values before sparse Jacobian preparation" begin
    checked_f!(du, u, p, t) = (
        any(x -> x <= 0, u) && throw(DomainError(u));
        du .= -u;
        nothing
    )
    alg = Rodas5P(
        autodiff = AutoSparse(
            AutoForwardDiff();
            sparsity_detector = TracerLocalSparsityDetector(),
            coloring_algorithm = GreedyColoringAlgorithm()
        )
    )
    integ = init(ODEProblem(checked_f!, [1.0, 1.0, 1.0], (0.0, 1.0)), alg)
    integ.u[3] = -99.0
    deleteat!(integ, 3)
    addat!(integ, 3, 2.0)
    @test integ.u == [1.0, 1.0, 2.0]
    @test all(v -> v[3] == 2.0, step_history(integ))
end

@testset "default cache hooks tolerate uninitialized sub-caches" begin
    integ = init(ODEProblem(decay!, [1.0, 1.0], (0.0, 1.0), [1.0, 1.0]))
    @test resize_non_user_cache!(integ, 2) === nothing
    @test deleteat_non_user_cache!(integ, 1:1) === nothing
    @test addat_non_user_cache!(integ, 3:3) === nothing
    @test integ.u == [1.0, 1.0]
end

function default_fbdf_integrator(n)
    p = fill(0.5, n + 2)
    p[(end - 1):end] .= 2000.0
    integ = init(
        ODEProblem(decay_prefix!, ones(n), (0.0, 0.002), p),
        DefaultODEAlgorithm(stiffalgfirst = true);
        abstol = 1.0e-8, reltol = 1.0e-8
    )
    full_cache(integ)
    fbdf_caches = filter(
        c -> c isa FBDFCache,
        default_subcaches(integ.cache)
    )
    @test length(fbdf_caches) == 2
    return integ, fbdf_caches
end

@testset "default solver resizes both FBDF nonlinear solvers" begin
    for (name, n, target, change!) in (
            ("addat!", 52, 54, integ -> addat!(integ, 53:54, [1.0, 1.0])),
            ("resize!", 52, 54, integ -> (resize!(integ, 54); integ.u[53:54] .= 1.0)),
            ("deleteat!", 54, 52, integ -> deleteat!(integ, 53:54)),
        )
        integ, fbdf_caches = default_fbdf_integrator(n)
        change!(integ)
        @test length(integ.u) == target
        @test all(fbdf_caches) do cache
            nlsolver = cache.nlsolver
            length(nlsolver.z) == length(nlsolver.tmp) == length(nlsolver.ztmp) == target &&
                length(nlsolver.cache.dz) == target &&
                size(nlsolver.cache.J) == (target, target) &&
                size(nlsolver.cache.W) == (target, target) &&
                length(nlsolver.cache.linsolve.b) == target
        end
        step!(integ)
    end
end

@testset "default solver callback $name resizes FBDF caches" for (name, n, target, change!) in (
        ("addat!", 52, 54, integ -> addat!(integ, 53:54, [1.0, 1.0])),
        ("resize!", 52, 54, integ -> (resize!(integ, 54); integ.u[53:54] .= 1.0)),
        ("deleteat!", 54, 52, integ -> deleteat!(integ, 53:54)),
    )
    p = fill(0.5, n + 2)
    p[(end - 1):end] .= 2000.0
    prob = ODEProblem(decay_prefix!, ones(n), (0.0, 0.002), p)
    cb = DiscreteCallback(
        (u, t, integ) -> t >= 0.001 && length(u) == n,
        integ -> change!(integ);
        save_positions = (false, false)
    )
    sol = solve(
        prob, DefaultODEAlgorithm(stiffalgfirst = true);
        callback = cb, tstops = [0.001], abstol = 1.0e-8, reltol = 1.0e-8
    )
    @test sol.retcode == ReturnCode.Success
    @test length(sol.u[end]) == target
end
