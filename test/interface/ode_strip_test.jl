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
vern_sol = solve(prob,Vern7())
@testset "Interpolation Stripping" begin
    @test isnothing(SciMLBase.strip_interpolation(rosenbrock_sol.interp).f)
    @test isnothing(SciMLBase.strip_interpolation(rosenbrock_sol.interp).cache.jac_config)
    @test isnothing(SciMLBase.strip_interpolation(rosenbrock_sol.interp).cache.grad_config)
end

@testset "Rosenbrock Solution Stripping" begin
    @test isnothing(SciMLBase.strip_solution(rosenbrock_sol).prob)
    @test isnothing(SciMLBase.strip_solution(rosenbrock_sol).alg)
    @test isnothing(SciMLBase.strip_solution(rosenbrock_sol).interp.f)
    @test isnothing(SciMLBase.strip_solution(rosenbrock_sol).interp.cache.jac_config)
    @test isnothing(SciMLBase.strip_solution(rosenbrock_sol).interp.cache.grad_config)
end 

@testset "TRBDF Solution Stripping" begin 
    @test isnothing(SciMLBase.strip_solution(TRBDF_sol).prob)
    @test isnothing(SciMLBase.strip_solution(TRBDF_sol).alg)
    @test isnothing(SciMLBase.strip_solution(TRBDF_sol).interp.f)
    @test isnothing(SciMLBase.strip_solution(TRBDF_sol).interp.cache.jac_config)
    @test isnothing(SciMLBase.strip_solution(TRBDF_sol).interp.cache.grad_config)
end


@test_throws SciMLBase.LazyInterpolationException SciMLBase.strip_solution(vern_sol)