using Test
using DiffEqBase
using FlexUnits
using FlexUnits.UnitRegistry
using DiffEqBase: ODE_DEFAULT_NORM

@testset "FlexUnits ODE_DEFAULT_NORM is RMS" begin
    # Issue #4610: array norms must match RMS on stripped values,
    # sqrt(mean(abs2(value(x)))), not sqrt(mean(abs(value(x)))).
    plain = [3.0, 4.0]
    expected = sqrt((abs2(3.0) + abs2(4.0)) / 2)
    @test ODE_DEFAULT_NORM(plain, 0.0) ≈ expected
    @test expected ≈ 3.5355339059327378

    u_arr = [3.0u"m", 4.0u"m"]
    @test ODE_DEFAULT_NORM(u_arr, 0.0) ≈ ODE_DEFAULT_NORM(plain, 0.0)
    @test ODE_DEFAULT_NORM(u_arr, 0.0) ≈ expected

    u_view = @view u_arr[:]
    @test u_view isa AbstractArray{<:FlexUnits.Quantity}
    @test !(u_view isa Array)
    @test ODE_DEFAULT_NORM(u_view, 0.0) ≈ ODE_DEFAULT_NORM(plain, 0.0)
    @test ODE_DEFAULT_NORM(u_view, 0.0) ≈ expected

    @test ODE_DEFAULT_NORM(3.0u"m", 0.0) == 3.0
end
