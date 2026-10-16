using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, Test
using OrdinaryDiffEqCore: AbstractControllerCache, IController, PIController,
    PIDController, PredictiveController, restore_controller!, snapshot_controller,
    sync_controllers!

scratch_fields(cache) = filter(n -> n ∉ (:controller, :EEst), fieldnames(typeof(cache)))

scratch(cache) = map(scratch_fields(cache)) do n
    v = getfield(cache, n)
    v isa AbstractArray ? copy(v) : v
end

function perturb!(cache)
    for n in scratch_fields(cache)
        v = getfield(cache, n)
        v isa AbstractArray ? (v .= 2 .* v .+ 1) : setfield!(cache, n, 2 * v + 1)
    end
    return cache
end

roundtrip_allocs(cache) = @allocated restore_controller!(cache, snapshot_controller(cache))

f!(du, u, p, t) = (du .= -u; nothing)
prob = ODEProblem(f!, [1.0, 2.0], (0.0, 1.0))

@testset "built-in $(nameof(typeof(controller))) cache" for controller in (
        IController(), PIController(7 // 50, 2 // 25), PIDController(0.7, -0.4),
        PredictiveController(),
    )
    cache = init(prob, Tsit5(); controller).controller_cache
    before = scratch(perturb!(cache))
    snap = snapshot_controller(cache)
    @test snap isa NamedTuple
    perturb!(cache)
    @test scratch(cache) != before
    @test restore_controller!(cache, snap) === nothing
    @test scratch(cache) == before

    # The snapshot must cover everything `sync_controllers!` transfers.
    other = perturb!(deepcopy(cache))
    sync_controllers!(other, cache)
    @test scratch(other) == scratch(cache)

    roundtrip_allocs(cache)
    @test roundtrip_allocs(cache) == 0
end

mutable struct HistoryControllerCache{E} <: AbstractControllerCache
    EEst::E
    history::Vector{Float64}
end

function OrdinaryDiffEqCore.sync_controllers!(
        dst::HistoryControllerCache, src::HistoryControllerCache
    )
    copyto!(dst.history, src.history)
    return nothing
end

@testset "custom cache falls back to deepcopy + sync_controllers!" begin
    cache = HistoryControllerCache(1.0, [1.0, 2.0, 3.0, 4.0])
    snap = snapshot_controller(cache)
    @test snap isa HistoryControllerCache
    @test snap.history !== cache.history
    cache.history .= 0
    @test restore_controller!(cache, snap) === nothing
    @test cache.history == [1.0, 2.0, 3.0, 4.0]
end
