using Test
using DiffEqBase
using FlexUnits
using FlexUnits.UnitRegistry
using Measurements
using MonteCarloMeasurements
using ReverseDiff
using Tracker
using DiffEqBase: ODE_DEFAULT_NORM

@testset "Extension ODE_DEFAULT_NORM is RMS" begin
    plain = [3.0, 4.0]
    expected = sqrt((abs2(3.0) + abs2(4.0)) / 2)
    @test ODE_DEFAULT_NORM(plain, 0.0) ≈ expected

    # Varied particles whose means are 3 and 4 (not constant fills).
    mcm3 = Particles([1.0, 2.0, 3.0, 4.0, 5.0])
    mcm4 = Particles([0.0, 2.0, 4.0, 6.0, 8.0])
    @test DiffEqBase.value(mcm3) == 3.0
    @test DiffEqBase.value(mcm4) == 4.0

    rd_tracked = ReverseDiff.track(plain)
    @test ODE_DEFAULT_NORM(rd_tracked, 0.0) ≈ expected  # TrackedArray positive control

    tr_tracked = Tracker.param(plain)
    @test ODE_DEFAULT_NORM(tr_tracked, 0.0) ≈ expected  # TrackedArray positive control

    cases = (
        (
            "FlexUnits",
            () -> [3.0u"m", 4.0u"m"],
            FlexUnits.Quantity,
            3.0u"m",
            3.0,
            nothing,
        ),
        (
            "Measurements",
            () -> [measurement(3.0, 0.0), measurement(4.0, 0.0)],
            Measurement,
            measurement(3.0, 0.0),
            3.0,
            nothing,
        ),
        (
            "MonteCarloMeasurements",
            () -> [mcm3, mcm4],
            AbstractParticles,
            mcm3,
            3.0,
            nothing,
        ),
        (
            "ReverseDiff",
            () -> [rd_tracked[1], rd_tracked[2]],
            ReverseDiff.TrackedReal,
            rd_tracked[1],
            3.0,
            ReverseDiff.track(0.0),
        ),
        (
            "Tracker",
            () -> [Tracker.TrackedReal(3.0), Tracker.TrackedReal(4.0)],
            Tracker.TrackedReal,
            Tracker.TrackedReal(3.0),
            3.0,
            Tracker.TrackedReal(0.0),
        ),
    )

    for (label, make_arr, ElT, scalar, scalar_expected, t_tracked) in cases
        @testset "$label" begin
            u_arr = make_arr()
            @test ODE_DEFAULT_NORM(u_arr, 0.0) ≈ ODE_DEFAULT_NORM(plain, 0.0)

            u_view = @view u_arr[:]
            @test u_view isa AbstractArray{<:ElT}
            @test !(u_view isa Array)
            @test ODE_DEFAULT_NORM(u_view, 0.0) ≈ ODE_DEFAULT_NORM(plain, 0.0)

            if t_tracked !== nothing
                @test ODE_DEFAULT_NORM(u_arr, t_tracked) ≈ expected
                @test ODE_DEFAULT_NORM(u_view, t_tracked) ≈ expected
            end

            @test ODE_DEFAULT_NORM(scalar, 0.0) == scalar_expected
        end
    end
end
