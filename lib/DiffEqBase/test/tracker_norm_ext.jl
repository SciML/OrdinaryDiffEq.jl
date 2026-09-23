using Test
using DiffEqBase
using Tracker
using DiffEqBase: ODE_DEFAULT_NORM

@testset "Tracker ODE_DEFAULT_NORM is RMS" begin
    # Issue #4610: AbstractArray{<:TrackedReal} norms must match RMS,
    # sqrt(mean(abs2(value(x)))), not sqrt(mean(abs(value(x)))).
    # TrackedArray methods already use abs2 and are a positive control.
    plain = [3.0, 4.0]
    expected = sqrt((abs2(3.0) + abs2(4.0)) / 2)
    @test ODE_DEFAULT_NORM(plain, 0.0) ≈ expected
    @test expected ≈ 3.5355339059327378

    tracked = Tracker.param(plain)
    @test ODE_DEFAULT_NORM(tracked, 0.0) ≈ expected

    u_arr = [Tracker.TrackedReal(3.0), Tracker.TrackedReal(4.0)]
    @test u_arr isa Array{<:Tracker.TrackedReal}
    @test ODE_DEFAULT_NORM(u_arr, 0.0) ≈ ODE_DEFAULT_NORM(plain, 0.0)
    @test ODE_DEFAULT_NORM(u_arr, 0.0) ≈ expected

    u_view = @view u_arr[:]
    @test u_view isa AbstractArray{<:Tracker.TrackedReal}
    @test !(u_view isa Array)
    @test ODE_DEFAULT_NORM(u_view, 0.0) ≈ expected

    t_tracked = Tracker.TrackedReal(0.0)
    @test ODE_DEFAULT_NORM(u_arr, t_tracked) ≈ expected
    @test ODE_DEFAULT_NORM(u_view, t_tracked) ≈ expected

    @test ODE_DEFAULT_NORM(Tracker.TrackedReal(3.0), 0.0) == 3.0
end
