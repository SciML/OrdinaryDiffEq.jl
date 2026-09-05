using OrdinaryDiffEqCore: IController, PIController
using OrdinaryDiffEq
using Reactant
using SciMLBase
using Test

f(u, p, t) = p .* u

function f!(du, u, p, t)
    du .= p .* u
    return nothing
end

struct CompiledODESolve{F, A, K}
    f::F
    alg::A
    kwargs::K
end

function (s::CompiledODESolve)(u, p)
    prob = ODEProblem(s.f, u, (0.0f0, 1.0f0), p)
    return solve(prob, s.alg; s.kwargs...)
end

u0 = Reactant.to_rarray(Float32[1, 2])
p0 = Reactant.to_rarray(Float32[-1])

@testset "Rejected first step" for rhs in (f, f!)
    maxiters_solver = CompiledODESolve(
        rhs, Tsit5(), (; maxiters = 1, dt = 1.0f0, abstol = eps(Float32), reltol = eps(Float32))
    )
    maxiters_sol = Reactant.@jit maxiters_solver(u0, p0)
    @test maxiters_sol.retcode == ReturnCode.MaxIters
    @test !SciMLBase.successful_retcode(maxiters_sol)
    @test maxiters_sol.t == Float32[0]
    @test Array(maxiters_sol.u[end]) == Float32[1, 2]
end

implicit_solver = CompiledODESolve(f, Rosenbrock23(), (;))
@test_throws ArgumentError Reactant.@jit implicit_solver(u0, p0)
fixed_solver_without_dt = CompiledODESolve(f, Tsit5(), (; adaptive = false))
@test_throws ArgumentError Reactant.@jit fixed_solver_without_dt(u0, p0)

solver_cases = (
    ("adaptive Tsit5", CompiledODESolve(f, Tsit5(), (;))),
    ("fixed Tsit5", CompiledODESolve(f, Tsit5(), (; adaptive = false, dt = 0.1f0))),
    ("in-place adaptive Tsit5", CompiledODESolve(f!, Tsit5(), (;))),
    ("custom PIController Tsit5", CompiledODESolve(f, Tsit5(), (; controller = PIController(0.14, 0.08)))),
    ("IController Tsit5", CompiledODESolve(f, Tsit5(), (; controller = IController(), abstol = 1.0f-7, reltol = 1.0f-5))),
    ("in-place fixed Tsit5", CompiledODESolve(f!, Tsit5(), (; adaptive = false, dt = 0.1f0))),
    ("adaptive Vern7", CompiledODESolve(f, Vern7(), (;))),
)

@testset "$name" for (name, solver) in solver_cases
    compiled = Reactant.compile(solver, (u0, p0))
    for rate in (0.0f0, -1.0f0, -3.0f0)
        sol = compiled(
            Reactant.to_rarray(Float32[1, 2]), Reactant.to_rarray(Float32[rate])
        )
        @test sol.u[end] isa Reactant.ConcreteRArray
        @test Array(sol.u[end]) ≈ Float32[exp(rate), 2exp(rate)] rtol = 5.0f-4
        @test sol.t == Float32[1]
        @test sol.retcode == ReturnCode.Success
        @test SciMLBase.successful_retcode(sol)
        @test sol.prob === nothing
        @test sol.stats === nothing
        @test sol.interp === nothing
    end
end
