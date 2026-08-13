using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqCore: DIRK, COEFFICIENT_MULTISTEP
using SciMLBase: DAEFunction
using Test

@testset "DAE _compute_rhs COEFFICIENT_MULTISTEP path (issue #4113)" begin
    f = DAEFunction{false}((du, u, p, t) -> du .- u)
    tmp = [0.5, -0.25]
    α = 0.75
    invγdt = 2.0
    uprev = [1.0, 2.0]
    z = [0.1, 0.2]
    tstep = 0.0
    p = nothing

    resid_dirk, u_dirk = OrdinaryDiffEqNonlinearSolve._compute_rhs(
        tmp, α, tstep, invγdt, p, uprev, f, z, DIRK
    )
    @test u_dirk ≈ uprev .+ z
    @test resid_dirk ≈ (tmp .+ α .* z) .* invγdt .- u_dirk

    resid_ms, u_ms = OrdinaryDiffEqNonlinearSolve._compute_rhs(
        tmp, α, tstep, invγdt, p, uprev, f, z, COEFFICIENT_MULTISTEP
    )
    @test u_ms ≈ z
    @test resid_ms ≈ (tmp .+ α .* z) .* invγdt .- z
    @test u_ms ≉ u_dirk
    @test resid_ms ≉ resid_dirk
end
