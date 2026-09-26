using Test
using DiffEqBase
using Measurements
using DiffEqBase: ODE_DEFAULT_NORM

@testset "Measurements ODE_DEFAULT_NORM is RMS" begin
    plain = [3.0, 4.0]
    expected = sqrt((abs2(3.0) + abs2(4.0)) / 2)
    @test ODE_DEFAULT_NORM(plain, 0.0) ≈ expected
    @test expected ≈ 3.5355339059327378

    u_arr = [measurement(3.0, 0.0), measurement(4.0, 0.0)]
    @test ODE_DEFAULT_NORM(u_arr, 0.0) ≈ ODE_DEFAULT_NORM(plain, 0.0)
    @test ODE_DEFAULT_NORM(u_arr, 0.0) ≈ expected

    u_view = @view u_arr[:]
    @test u_view isa AbstractArray{<:Measurement}
    @test !(u_view isa Array)
    @test ODE_DEFAULT_NORM(u_view, 0.0) ≈ ODE_DEFAULT_NORM(plain, 0.0)
    @test ODE_DEFAULT_NORM(u_view, 0.0) ≈ expected

    @test ODE_DEFAULT_NORM(measurement(3.0, 0.0), 0.0) == 3.0
end
