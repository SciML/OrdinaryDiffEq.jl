using OrdinaryDiffEq, Test
import SciMLBase

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 0.5)
prob = ODEProblem(lorenz!, u0, tspan)

rosenbrock_sol = solve(prob, Rosenbrock23())
TRBDF_sol = solve(prob, TRBDF2())
vern_sol = solve(prob, Vern7())
default_sol = solve(prob)
@testset "Interpolation Stripping" begin
    @test isnothing(SciMLBase.strip_interpolation(rosenbrock_sol.interp).f)
    @test isnothing(SciMLBase.strip_interpolation(rosenbrock_sol.interp).cache.jac_config)
    @test isnothing(SciMLBase.strip_interpolation(rosenbrock_sol.interp).cache.grad_config)
end

@testset "Rosenbrock Solution Stripping" begin
    stripped_sol = SciMLBase.strip_solution(rosenbrock_sol)
    @test stripped_sol.prob isa NamedTuple
    @test isnothing(SciMLBase.strip_solution(rosenbrock_sol, strip_alg = true).alg)
    @test isnothing(stripped_sol.interp.f)
    @test isnothing(stripped_sol.interp.cache.jac_config)
    @test isnothing(stripped_sol.interp.cache.grad_config)
    @test isnothing(stripped_sol.interp.cache.uf)
    @test isnothing(stripped_sol.interp.cache.tf)
end

@testset "TRBDF Solution Stripping" begin
    stripped_sol = SciMLBase.strip_solution(TRBDF_sol)
    @test stripped_sol.prob isa NamedTuple
    @test isnothing(SciMLBase.strip_solution(TRBDF_sol, strip_alg = true).alg)
    @test isnothing(stripped_sol.interp.f)
    @test isnothing(stripped_sol.interp.cache.nlsolver)
end

@testset "Default Solution Stripping" begin
    stripped_sol = SciMLBase.strip_solution(default_sol)
    @test isnothing(stripped_sol.interp.cache.args)
end

@test_throws SciMLBase.LazyInterpolationException SciMLBase.strip_solution(vern_sol)
