using DelayDiffEq, Test
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqVerner
using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqRosenbrock
using SciMLBase: ReturnCode
using ADTypes: AutoFiniteDiff
using DiffEqBase: full_cache

ddef!(du, u, h, p, t) = (du .= -50 .* (u .- 0.5 .* h(p, t - 0.1; idxs = 1)); nothing)
h0(p, t; idxs = nothing) = idxs === nothing ? ones(3) : 1.0
const prob = DDEProblem(ddef!, ones(3), h0, (0.0, 1.0); constant_lags = [0.1])

function state_buffers(integ)
    bufs = Any[integ.u, integ.uprev, integ.uprev2, integ.integrator.u]
    append!(bufs, integ.k, integ.integrator.k, full_cache(integ.cache))
    for c in (integ.fsalfirst, integ.fsallast)
        c isa AbstractArray && push!(bufs, c)
    end
    return bufs
end

@testset "deleteat!/addat! with $(nameof(typeof(alg)))" for alg in (
        Tsit5(), Vern7(), BS3(), Trapezoid(autodiff = AutoFiniteDiff()),
        KenCarp4(autodiff = AutoFiniteDiff()), Rosenbrock23(autodiff = AutoFiniteDiff()),
    )
    integ = init(prob, MethodOfSteps(alg))
    for _ in 1:3
        step!(integ)
    end
    u = copy(integ.u)
    deleteat!(integ, 2)
    @test integ.u == u[[1, 3]]
    @test all(c -> length(c) == 2, state_buffers(integ))
    solve!(integ)
    @test integ.sol.retcode == ReturnCode.Success

    integ = init(prob, MethodOfSteps(alg))
    for _ in 1:3
        step!(integ)
    end
    addat!(integ, 4:4)
    integ.u[4] = 1.0
    u_modified!(integ, true)
    @test all(c -> length(c) == 4, state_buffers(integ))
    solve!(integ)
    @test integ.sol.retcode == ReturnCode.Success
end

# A lag of 8 fixed steps puts every propagated discontinuity on the step grid of both
# integrators, so the continued and the fresh solve take identical steps and differ only
# by rounding (largest observed: 1.3e-13 for Rodas5P with a finite-difference Jacobian).
@testset "deleteat! matches a fresh solve of the reduced system with $(nameof(typeof(alg)))" for alg in (
        Tsit5(), Vern7(), BS3(), Trapezoid(autodiff = AutoFiniteDiff()),
        KenCarp4(autodiff = AutoFiniteDiff()), Rodas5P(autodiff = AutoFiniteDiff()),
    )
    u0 = [1.0, 2.0, 3.0]
    hist(p, t; idxs = nothing) = idxs === nothing ? u0 : u0[idxs]
    lag = 1 / 8
    f!(du, u, h, p, t) = (du .= -50 .* (u .- 0.5 .* h(p, t - lag; idxs = 1)); nothing)
    kwargs = (; adaptive = false, dt = 1 / 64)
    integ = init(
        DDEProblem(f!, u0, hist, (0.0, 1.0); constant_lags = [lag]), MethodOfSteps(alg);
        kwargs...
    )
    for _ in 1:4
        step!(integ)
    end
    kept = [1, 3]
    pre = deepcopy(integ.sol)
    href(p, t; idxs = nothing) = t < 0 ? hist(p, t; idxs = kept[something(idxs, :)]) :
        idxs === nothing ? pre(t)[kept] : pre(t; idxs = kept[idxs])
    ref = solve(
        DDEProblem(f!, integ.u[kept], href, (integ.t, 1.0); constant_lags = [lag]),
        MethodOfSteps(alg); kwargs...
    )
    deleteat!(integ, 2)
    solve!(integ)
    @test integ.sol.retcode == ReturnCode.Success
    @test integ.t == ref.t[end] == 1.0
    @test integ.u ≈ ref.u[end] rtol = 1.0e-12
end

@testset "grown step history is zeroed" begin
    # Leftover array capacity stands in for uninitialized memory; this relies on
    # shrinking a Vector keeping its contents.
    poison = 2.4e77
    v = resize!(fill!(Vector{Float64}(undef, 5), poison), 3)
    resize!(v, 5)
    @test all(==(poison), v[4:5])

    history(integ) = (integ.uprev, integ.uprev2, integ.cache.uprev3)
    shrink_grow = (
        "resize!" => integ -> (resize!(integ, 3); resize!(integ, 5)),
        "deleteat!/addat!" => integ -> (deleteat!(integ, 4:5); addat!(integ, 4:5)),
    )
    @testset "$name" for (name, shrink_grow!) in shrink_grow
        integ = init(prob, MethodOfSteps(Trapezoid()))
        for _ in 1:4
            step!(integ)
        end
        resize!(integ, 5)
        foreach(v -> v[4:5] .= poison, history(integ))
        shrink_grow!(integ)
        @test all(v -> length(v) == 5 && iszero(v[4:5]), history(integ))
        integ.u[4:5] .= 1.0
        u_modified!(integ, true)
        solve!(integ)
        @test integ.sol.retcode == ReturnCode.Success
    end
end
