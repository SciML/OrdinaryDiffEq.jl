using OrdinaryDiffEqCore, Test

function popall!(queue)
    values = eltype(queue)[]
    while !isempty(queue)
        push!(values, pop!(queue))
    end
    return values
end

@testset "Time queue initialization" begin
    @test popall!(
        OrdinaryDiffEqCore.initialize_tstops(
            Float64, (0.75, 0.25, -1.0, 2.0), (0.5,), (0.0, 1.0)
        )
    ) == [0.25, 0.5, 0.75, 1.0]

    @test popall!(
        OrdinaryDiffEqCore.initialize_tstops(
            Float64, (0.75, 0.25, -1.0, 2.0), (0.5,), (1.0, 0.0)
        )
    ) == [-0.75, -0.5, -0.25, 0.0]

    @test popall!(OrdinaryDiffEqCore.initialize_saveat(Float64, 0.25, (0.0, 1.0))) ==
        [0.25, 0.5, 0.75, 1.0]
    @test popall!(
        OrdinaryDiffEqCore.initialize_saveat(Float64, (0.25, 1.0, 2.0), (1.0, 0.0))
    ) == [-0.25]

    @test popall!(
        OrdinaryDiffEqCore.initialize_d_discontinuities(
            Float64, (-1.0, 0.0, 0.25, 1.0, 2.0), (0.0, 1.0)
        )
    ) == [0.0, 0.25, 1.0, 2.0]
    @test popall!(
        OrdinaryDiffEqCore.initialize_d_discontinuities(
            Float64, (2.0, 1.0, 0.75, 0.0, -1.0), (1.0, 0.0)
        )
    ) == [-1.0, -0.75, 0.0, 1.0]

    discontinuities = OrdinaryDiffEqCore.initialize_d_discontinuities(
        Float64, (-1.0,), (0.0, 1.0)
    )
    OrdinaryDiffEqCore.reinit_d_discontinuities!(
        Float64, discontinuities, (-1.0, 0.0, 2.0), (0.0, 1.0)
    )
    @test popall!(discontinuities) == [0.0, 2.0]

    float32_tstops = OrdinaryDiffEqCore.initialize_tstops(
        Float32, (Float32(0.5),), (), (Float32(0), Float32(1))
    )
    @test popall!(float32_tstops) == Float32[0.5, 1.0]
end
