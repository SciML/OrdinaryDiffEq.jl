using OrdinaryDiffEq, Test
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqBDF: QNDF
using OrdinaryDiffEqFIRK: RadauIIA5, AdaptiveRadau
using OrdinaryDiffEqLowOrderRK: BS5
import SciMLBase

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    return
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

@testset "Composite Solution Stripping" begin
    composite_sol = solve(prob, AutoTsit5(Rosenbrock23()))
    @test all(==(1), composite_sol.alg_choice)
    stripped_sol = SciMLBase.strip_solution(composite_sol)
    @test isnothing(stripped_sol.interp.f)
    @test isnothing(stripped_sol.interp.cache.choice_function)
    @test isnothing(stripped_sol.interp.cache.caches[2].jac_config)
    @test stripped_sol(0.33) == composite_sol(0.33)
    @test SciMLBase.strip_solution(solve(prob, AutoVern7(Rosenbrock23(); lazy = false))) isa ODESolution
end

@testset "Stripped solutions interpolate like the original" begin
    lorenz(u, p, t) = (du = similar(u); lorenz!(du, u, p, t); du)
    ts = range(tspan...; length = 101)
    algs = (
        Rosenbrock23(), Rodas5P(), FBDF(), QNDF(), RadauIIA5(), AdaptiveRadau(),
        BS5(lazy = false), DefaultODEAlgorithm(lazy = Val(false)),
        AutoTsit5(Rosenbrock23()), AutoVern7(Rodas5P(); lazy = false),
    )
    for f in (lorenz!, lorenz), alg in algs
        sol = solve(ODEProblem(f, u0, tspan), alg)
        stripped_sol = SciMLBase.strip_solution(sol)
        @test isnothing(stripped_sol.interp.f)
        @test all(isapprox(stripped_sol(t), sol(t); rtol = 1.0e-12) for t in ts)
    end

    holds(x, y, seen = IdDict()) = x === y || (
        !(x isa Union{Module, Type}) && !haskey(seen, x) && (seen[x] = true;
            x isa Array ? any(i -> isassigned(x, i) && holds(x[i], y, seen), eachindex(x)) :
            any(i -> isdefined(x, i) && holds(getfield(x, i), y, seen), 1:nfields(x)))
    )
    parameter_rate(t) = 1.0
    growth = ODEProblem((du, u, p, t) -> (du .= p.rate(t) .* u), [1.0], (0.0, 1.0), (rate = parameter_rate,))
    ts = range(0.013, 0.987; length = 50)
    for alg in (AdaptiveRadau(), RadauIIA5(), DefaultODEAlgorithm(lazy = Val(false)))
        sol = solve(growth, alg; reltol = 1.0e-10, abstol = 1.0e-10)
        stripped_sol = SciMLBase.strip_solution(sol; strip_alg = true)
        @test !holds(stripped_sol.interp, parameter_rate)
        @test all(stripped_sol(t) == sol(t) for t in ts)
    end
    sol = SciMLBase.strip_solution(solve(growth, DefaultODEAlgorithm(lazy = Val(false)); reltol = 1.0e-10, abstol = 1.0e-10))
    @test maximum(t -> abs(sol(t)[1] - exp(t)), ts) < 1.0e-11
end
