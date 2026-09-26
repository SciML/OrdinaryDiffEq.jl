using OrdinaryDiffEqCore
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqBDF
using Unitful
using Test

function expected_default_tolerances(u::AbstractArray{T}) where {T <: Number}
    abstol = real(convert(T, oneunit(T) * 1 // 10^6))
    reltol = real(convert(T, oneunit(T) * 1 // 10^3))
    return abstol, reltol
end

function decay!(du, u, p, t)
    du[1] = -u[1]
    return du[2] = -2u[2]
end

function decay(u, p, t)
    return [-u[1], -2u[2]]
end

@testset "default tolerance resolution preserves values and types" begin
    for (T, alg) in (
            (Float64, Tsit5()),
            (Float32, Tsit5()),
            (BigFloat, Tsit5()),
            (Float64, Rosenbrock23()),
        )
        u0 = T[1, 2]
        tspan = (zero(T), one(T))
        prob = ODEProblem(decay!, u0, tspan)
        exp_abstol, exp_reltol = expected_default_tolerances(u0)

        resolved = OrdinaryDiffEqCore.resolve_ode_tolerances(prob, u0, nothing, nothing)
        @test resolved[1] == exp_abstol
        @test resolved[2] == exp_reltol
        @test typeof(resolved[1]) == typeof(exp_abstol)
        @test typeof(resolved[2]) == typeof(exp_reltol)

        integ_default = init(prob, alg)
        integ_explicit = init(prob, alg; abstol = exp_abstol, reltol = exp_reltol)
        @test integ_default.opts.abstol == integ_explicit.opts.abstol
        @test integ_default.opts.reltol == integ_explicit.opts.reltol
        @test typeof(integ_default.opts.abstol) == typeof(integ_explicit.opts.abstol)
        @test typeof(integ_default.opts.reltol) == typeof(integ_explicit.opts.reltol)

        sol_default = solve(prob, alg)
        sol_explicit = solve(prob, alg; abstol = exp_abstol, reltol = exp_reltol)
        @test sol_default.retcode == sol_explicit.retcode
        @test sol_default.stats.naccept == sol_explicit.stats.naccept
        @test sol_default.stats.nreject == sol_explicit.stats.nreject
        @test sol_default.t == sol_explicit.t
        @test sol_default.u == sol_explicit.u
    end
end

@testset "vector abstol still works" begin
    u0 = [1.0, 2.0]
    prob = ODEProblem(decay!, u0, (0.0, 1.0))
    abstol = [1.0e-8, 1.0e-6]
    reltol = [1.0e-6, 1.0e-4]
    sol = solve(prob, Tsit5(); abstol, reltol)
    @test SciMLBase.successful_retcode(sol)
    integ = init(prob, Tsit5(); abstol, reltol)
    @test integ.opts.abstol == abstol
    @test integ.opts.reltol == reltol
end

@testset "oop Float64 default equality" begin
    prob = ODEProblem(decay, [1.0, 2.0], (0.0, 1.0))
    exp_abstol, exp_reltol = expected_default_tolerances(prob.u0)
    sol_a = solve(prob, Tsit5())
    sol_b = solve(prob, Tsit5(); abstol = exp_abstol, reltol = exp_reltol)
    @test sol_a.stats.naccept == sol_b.stats.naccept
    @test sol_a.u == sol_b.u
end

@testset "resolve leaves user Float64 tolerances unchanged in type" begin
    u0 = [1.0, 2.0]
    prob = ODEProblem(decay!, u0, (0.0, 1.0))
    a, r = OrdinaryDiffEqCore.resolve_ode_tolerances(prob, u0, 1.0e-6, 1.0e-6)
    @test a === 1.0e-6
    @test r === 1.0e-6
    @test typeof(a) == Float64
    a0, r0 = OrdinaryDiffEqCore.resolve_ode_tolerances(prob, u0, nothing, nothing)
    @test typeof(a0) == typeof(a)
    @test typeof(r0) == typeof(r)
end

@testset "problem kwargs abstol takes precedence when call omits it" begin
    u0 = [1.0, 2.0]
    prob = ODEProblem(decay!, u0, (0.0, 1.0); abstol = 1.0e-8, reltol = 1.0e-5)
    integ = init(prob, Tsit5())
    @test integ.opts.abstol == 1.0e-8
    @test integ.opts.reltol == 1.0e-5
    sol_a = solve(prob, Tsit5())
    sol_b = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-5)
    @test sol_a.stats.naccept == sol_b.stats.naccept
    @test sol_a.u == sol_b.u
end

@testset "unitful u0 default tolerances" begin
    u0 = [1.0u"m", 2.0u"m"]
    tspan = (0.0u"s", 1.0u"s")
    function decay_unitful!(du, u, p, t)
        du[1] = -u[1] / 1u"s"
        return du[2] = -2u[2] / 1u"s"
    end
    prob = ODEProblem(decay_unitful!, u0, tspan)
    # Master defaults for this problem (vector Quantity{Float64} in meters).
    exp_a = [1.0e-6u"m", 1.0e-6u"m"]
    exp_r = [1.0e-3u"m", 1.0e-3u"m"]
    integ = init(prob, Tsit5())
    @test integ.opts.abstol == exp_a
    @test integ.opts.reltol == exp_r
    @test typeof(integ.opts.abstol) == typeof(exp_a)
    @test typeof(integ.opts.reltol) == typeof(exp_r)
    sol_a = solve(prob, Tsit5())
    sol_b = solve(prob, Tsit5(); abstol = exp_a, reltol = exp_r)
    @test sol_a.stats.naccept == sol_b.stats.naccept
    @test sol_a.u == sol_b.u
end

@testset "Float32 DAE default tolerances (DFBDF)" begin
    # index-1 semi-explicit DAE: u' = -u, 0 = u - v
    function f!(out, du, u, p, t)
        out[1] = du[1] + u[1]
        return out[2] = u[1] - u[2]
    end
    u0 = Float32[1, 1]
    du0 = Float32[-1, 0]
    tspan = (0.0f0, 1.0f0)
    prob = DAEProblem(f!, du0, u0, tspan; differential_vars = [true, false])
    exp_a, exp_r = expected_default_tolerances(u0)
    resolved = OrdinaryDiffEqCore.resolve_ode_tolerances(prob, u0, nothing, nothing)
    @test resolved[1] == exp_a
    @test resolved[2] == exp_r
    integ = init(prob, DFBDF())
    @test integ.opts.abstol == exp_a
    @test integ.opts.reltol == exp_r
    sol_a = solve(prob, DFBDF())
    sol_b = solve(prob, DFBDF(); abstol = exp_a, reltol = exp_r)
    @test sol_a.stats.naccept == sol_b.stats.naccept
    @test sol_a.u == sol_b.u
end

@testset "SDE-resolved concrete tolerances pass through unchanged" begin
    # StochasticDiffEq resolves to 1//10^2 before _ode_init; Core must not rewrite them.
    u0 = [1.0, 2.0]
    prob = ODEProblem(decay!, u0, (0.0, 1.0)) # type only matters for Discrete check
    sde_abstol = 1.0e-2
    sde_reltol = 1.0e-2
    a, r = OrdinaryDiffEqCore.resolve_ode_tolerances(prob, u0, sde_abstol, sde_reltol)
    @test a === sde_abstol
    @test r === sde_reltol
end

function count_kwbody_specializations(mod, name::Symbol)
    n = 0
    prefix = "#" * string(name) * "#"
    for s in names(mod; all = true)
        startswith(string(s), prefix) || continue
        fn = getfield(mod, s)
        fn isa Function || continue
        for m in methods(fn).ms
            n += count(!isnothing, Base.specializations(m))
        end
    end
    return n
end

@testset "tol solve does not add _ode_init kwbody specialization" begin
    # Property this PR exists for: after a default solve warms the Float64 path,
    # an explicit abstol/reltol solve must not compile a new heavy `_ode_init`
    # keyword body (SciMLBase.__init may grow; the heavy body must not).
    function f!(du, u, p, t)
        du .= -u
        return nothing
    end
    prob = ODEProblem(f!, [1.0, 2.0], (0.0, 1.0))
    solve(prob, Tsit5())
    n0 = count_kwbody_specializations(OrdinaryDiffEqCore, :_ode_init)
    solve(prob, Tsit5(); abstol = 1.0e-6, reltol = 1.0e-6)
    n1 = count_kwbody_specializations(OrdinaryDiffEqCore, :_ode_init)
    @test n1 == n0
    @test n0 > 0
end

@testset "ODE solve inference with and without tolerances" begin
    prob = ODEProblem(decay, [1.0, 2.0], (0.0, 1.0))
    @inferred solve(prob, Tsit5())
    @inferred solve(prob, Tsit5(); abstol = 1.0e-6, reltol = 1.0e-6)
end
