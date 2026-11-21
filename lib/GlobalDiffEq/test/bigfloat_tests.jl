using GlobalDiffEq, OrdinaryDiffEqTsit5
using OrdinaryDiffEqSSPRK
using Test

@testset "BigFloat support" begin
    function f_bf!(du, u, p, t)
        du[1] = -u[1]
    end

    u0_bf = BigFloat[1.0]
    tspan_bf = (BigFloat(0.0), BigFloat(1.0))
    prob_bf = ODEProblem(f_bf!, u0_bf, tspan_bf)

    sol_bf = solve(
        prob_bf, GlobalRichardson(SSPRK33()),
        dt = BigFloat(0.1), reltol = BigFloat(1.0e-3), abstol = BigFloat(1.0e-6)
    )

    @test eltype(sol_bf.u[end]) == BigFloat
    # Check solution is reasonable (e^-1 ≈ 0.368)
    @test isapprox(sol_bf.u[end][1], exp(BigFloat(-1)), rtol = 1.0e-4)
end
