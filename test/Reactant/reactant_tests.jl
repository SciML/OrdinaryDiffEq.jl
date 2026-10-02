using OrdinaryDiffEqCore: IController, PIController, PIDController
using OrdinaryDiffEq
using OrdinaryDiffEqCore
using DiffEqBase
using FastPower
using Logging
using Reactant
using Random
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
    ("PIDController Tsit5", CompiledODESolve(f, Tsit5(), (; controller = PIDController(0.7, -0.4), abstol = 1.0f-7, reltol = 1.0f-5))),
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

@testset "Fixed-step endpoint clipping" for rhs in (f, f!), direction in (1.0f0, -1.0f0)
    function step_count(u, p)
        integrator = init(
            ODEProblem(rhs, u, (0.0f0, direction), p), Tsit5();
            adaptive = false, dt = direction * 0.01f0, save_everystep = false
        )
        solve!(integrator)
        return integrator.iter
    end
    @test Reactant.@jit(step_count(u0, p0)) == step_count(Float32[1, 2], Float32[-1])
end

@testset "Initial-step NaN fallback" for rhs in (f, f!)
    function first_dt(u, p)
        return init(ODEProblem(rhs, u, (0.0f0, 1.0f0), p), Tsit5(); dtmin = 1.0f-5).dt
    end
    nan_p = Reactant.to_rarray(Float32[NaN])
    @test Reactant.@jit(first_dt(u0, nan_p)) == first_dt(Float32[1, 2], Float32[NaN])
end

unstable_check_solver = CompiledODESolve(
    f, Tsit5(), (; adaptive = false, dt = 0.1f0, unstable_check = (dt, u, p, t) -> false)
)
@test_throws ArgumentError Reactant.@jit unstable_check_solver(u0, p0)

@testset "Compiled failure retcodes match host" for (name, kwargs) in (
        ("fixed Tsit5", (; adaptive = false, dt = 0.1f0)),
        ("adaptive Tsit5", (;)),
    )
    solver = CompiledODESolve(f, Tsit5(), kwargs)
    u = Float32[1]
    compiled = Reactant.compile(
        solver, (Reactant.to_rarray(u), Reactant.to_rarray(Float32[-1]))
    )
    for p in (Float32[-1], Float32[NaN], Float32[1.0f30])
        host = solver(u, p)
        traced = compiled(Reactant.to_rarray(u), Reactant.to_rarray(p))
        @test traced.retcode == host.retcode
        if p == Float32[-1]
            @test host.retcode == ReturnCode.Success
            @test Array(traced.u[end]) ≈ Float32[0.3678795] rtol = 5.0f-4
        else
            @test !SciMLBase.successful_retcode(host)
            @test !Bool(SciMLBase.successful_retcode(traced))
        end
    end
end

@testset "Compiled failure retcodes equal the host's" begin
    u = Float32[1]
    cases = (
        ("small first step", (1.5f0, 2.0f0), (; dt = 1.5f0 * eps(1.5f0)), ((Float32[-1], ReturnCode.Success),)),
        (
            "dtmin after rejected steps", (0.0f0, 1.0f0), (; dt = 0.1f0, dtmin = 1.0f-5),
            ((Float32[NaN], ReturnCode.DtLessThanMin), (Float32[1.0f30], ReturnCode.DtLessThanMin), (Float32[-1], ReturnCode.Success)),
        ),
        ("NaN dt", (0.0f0, 1.0f0), (; dt = NaN32), ((Float32[-1], ReturnCode.DtNaN),)),
        (
            "PIDController", (0.0f0, 1.0f0), (; controller = PIDController(0.7, -0.4), maxiters = 1000),
            ((Float32[-1], ReturnCode.Success), (Float32[NaN], ReturnCode.DtNaN), (Float32[1.0f30], ReturnCode.DtNaN)),
        ),
    )
    @testset "$name" for (name, tspan, kwargs, params) in cases
        solver = (u, p) -> solve(ODEProblem(f, u, tspan, p), Tsit5(); kwargs...)
        compiled = Reactant.compile(
            solver, (Reactant.to_rarray(u), Reactant.to_rarray(Float32[-1]))
        )
        for (p, expected) in params
            host = solver(u, p)
            traced = compiled(Reactant.to_rarray(u), Reactant.to_rarray(p))
            @test host.retcode == expected
            @test traced.retcode == host.retcode
        end
    end

    @testset "existing failure code, adaptive = $adaptive" for adaptive in (false, true)
        function solve_failed(u, p)
            integrator = init(
                ODEProblem(f, u, (0.0f0, 1.0f0), p), Tsit5(); adaptive, dt = 0.1f0
            )
            integrator.sol = SciMLBase.solution_new_retcode(integrator.sol, ReturnCode.Unstable)
            return solve!(integrator)
        end
        host = solve_failed(u, Float32[-1])
        traced = Reactant.@jit solve_failed(Reactant.to_rarray(u), Reactant.to_rarray(Float32[-1]))
        @test host.retcode == ReturnCode.Unstable
        @test traced.retcode == host.retcode
    end
end

@testset "Staged error check equals the host check" begin
    prob = ODEProblem(f, Float32[1], (0.0f0, 1.0f0), Float32[-1])
    nchecked = 0
    for adaptive in (true, false), dt in (NaN32, 1.0f-30, 1.0f-6, 0.1f0),
            iter in (1, 11), accept_step in (false, true), unew in (1.0f0, NaN32),
            last_stepfail in (false, true), stored in (ReturnCode.Default, ReturnCode.MaxIters)

        integrator = init(prob, Tsit5(); adaptive, dt = 0.1f0, dtmin = 1.0f-5, maxiters = 10)
        integrator.dt = dt
        integrator.iter = iter
        integrator.accept_step = accept_step
        integrator.u .= unew
        integrator.last_stepfail = last_stepfail
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol, stored)
        host = with_logger(() -> SciMLBase.check_error(integrator), NullLogger())
        failed, code = DiffEqBase.staged_check_error(integrator)
        @test code == host
        @test failed == (host != ReturnCode.Success)
        nchecked += 1
    end
    @test nchecked == 256
end

@testset "Compiled eps(t) is exact" begin
    Random.seed!(0)
    below(dt, t) = OrdinaryDiffEqCore.dt_below_time_eps.(dt, t)
    @testset "$T" for T in (Float32, Float64)
        U = Base.uinttype(T)
        ts = T[]
        for biased in 0:(Base.exponent_mask(T) >> Base.significand_bits(T))
            x = reinterpret(T, U(biased) << Base.significand_bits(T))
            append!(ts, (x, nextfloat(x), prevfloat(max(x, nextfloat(zero(T))))))
        end
        append!(ts, (zero(T), floatmin(T), prevfloat(floatmin(T)), floatmax(T), T(NaN)))
        append!(ts, reinterpret.(T, rand(U, 2000)))
        append!(ts, -ts)
        epss = eps.(ts)
        for dts in (epss, nextfloat.(epss), prevfloat.(epss), -epss, zero(ts), reinterpret.(T, rand(U, length(ts))))
            traced = Array(Reactant.@jit below(Reactant.to_rarray(dts), Reactant.to_rarray(ts)))
            @test traced == (abs.(dts) .<= abs.(eps.(ts)))
        end
        # Compiled code flushes subnormals, so a subnormal spacing is raised to `floatmin`.
        spacing(t) = OrdinaryDiffEqCore.value_eps.(t)
        traced = Array(Reactant.@jit spacing(Reactant.to_rarray(ts)))
        @test isequal(traced, map(t -> isfinite(t) ? max(eps(t), floatmin(T)) : T(NaN), ts))
    end
end

@testset "Compiled controller power follows FastPower" begin
    fp(x, y) = OrdinaryDiffEqCore.controller_fastpower(x, y)
    @testset "$T" for T in (Float32, Float64)
        compiled = Reactant.compile(
            fp, (Reactant.ConcreteRNumber(one(T)), Reactant.ConcreteRNumber(one(T)))
        )
        Random.seed!(1)
        xs = T[0, NaN, Inf, 1.0e30, floatmax(T), floatmin(T), 1.0e-30, 0.5, 1, 1.5, 2, 100]
        ys = T[0.14, -0.08, 0.7, -0.4, 1 // 6, 0.2, Inf, 0]
        pairs = vcat(vec(collect(Iterators.product(xs, ys))), [(T(exp(20randn())), T(rand() - 0.5)) for _ in 1:200])
        for (x, y) in pairs
            host = FastPower.fastpower(x, y)
            traced = Float64(compiled(Reactant.ConcreteRNumber(x), Reactant.ConcreteRNumber(y)))
            # Both round an `exp2` evaluated in `Float32`, by different implementations.
            @test isfinite(traced) == isfinite(host)
            @test isapprox(traced, host; rtol = 2 * eps(Float32), nans = true)
        end
    end
end
