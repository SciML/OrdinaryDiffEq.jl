using StochasticDiffEqHighOrder
using SciMLBase
using Test

@testset "Float32 SOSRA OOP stages stay Float32" begin
    seen = DataType[]
    f(u, p, t) = (push!(seen, eltype(u)); -u)
    g(u, p, t) = (push!(seen, eltype(u)); 0.1f0 .+ zero(u))
    prob = SDEProblem{false}(f, g, Float32[1, 0.5, 0.25, 0.125], (0.0f0, 1.0f0))
    sol = solve(prob, SOSRA(); seed = 1)
    @test sol.retcode == ReturnCode.Success
    @test eltype(sol.u[end]) === Float32
    @test !isempty(seen)
    @test all(T -> T === Float32, seen)
end

@testset "SRI() OOP solves Vector u" begin
    f(u, p, t) = -u
    g(u, p, t) = 0.1 .* u
    prob = SDEProblem{false}(f, g, [1.0, 0.5], (0.0, 1.0))
    sol = solve(prob, SRI(); seed = 1)
    @test sol.retcode == ReturnCode.Success
    @test length(sol.u[end]) == 2
    @test all(isfinite, sol.u[end])
end
