using Test
using DiffEqBase
using Unitful
using DiffEqBase: ODE_DEFAULT_NORM

@testset "Unitful ODE_DEFAULT_NORM is RMS" begin
    # Issue #4596: Unitful array norms must match RMS on stripped values,
    # sqrt(mean(abs2(value(x)))), not sqrt(mean(abs(value(x)))).
    # Reference: RMS of [3, 4] is sqrt((9+16)/2) ≈ 3.5355339059327378
    plain = [3.0, 4.0]
    expected = sqrt((abs2(3.0) + abs2(4.0)) / 2)
    @test ODE_DEFAULT_NORM(plain, 0.0) ≈ expected
    @test expected ≈ 3.5355339059327378

    u_arr = [3.0, 4.0]u"m"
    @test ODE_DEFAULT_NORM(u_arr, 0.0) ≈ ODE_DEFAULT_NORM(plain, 0.0)
    @test ODE_DEFAULT_NORM(u_arr, 0.0) ≈ expected

    # SubArray hits the AbstractArray{<:Quantity} method, not Array
    u_view = @view u_arr[:]
    @test u_view isa AbstractArray{<:Unitful.AbstractQuantity}
    @test !(u_view isa Array)
    @test ODE_DEFAULT_NORM(u_view, 0.0) ≈ ODE_DEFAULT_NORM(plain, 0.0)
    @test ODE_DEFAULT_NORM(u_view, 0.0) ≈ expected

    @test ODE_DEFAULT_NORM(3.0u"m", 0.0) == 3.0
end
