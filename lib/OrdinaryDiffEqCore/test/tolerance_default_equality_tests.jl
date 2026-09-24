using OrdinaryDiffEqCore
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqRosenbrock
using DiffEqBase
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

        resolved = DiffEqBase.resolve_ode_tolerances(prob, u0, nothing, nothing)
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
    a, r = DiffEqBase.resolve_ode_tolerances(prob, u0, 1.0e-6, 1.0e-6)
    @test a === 1.0e-6
    @test r === 1.0e-6
    @test typeof(a) == Float64
    a0, r0 = DiffEqBase.resolve_ode_tolerances(prob, u0, nothing, nothing)
    @test typeof(a0) == typeof(a)
    @test typeof(r0) == typeof(r)
end

@testset "with_resolved_ode_tolerances unifies kwargs types" begin
    u0 = [1.0, 2.0]
    prob = ODEProblem(decay!, u0, (0.0, 1.0))
    kw_none = DiffEqBase.with_resolved_ode_tolerances(prob, u0, (;))
    kw_tol = DiffEqBase.with_resolved_ode_tolerances(
        prob, u0, (; abstol = 1.0e-6, reltol = 1.0e-6)
    )
    @test typeof(kw_none.abstol) == typeof(kw_tol.abstol)
    @test typeof(kw_none.reltol) == typeof(kw_tol.reltol)
end

@testset "problem kwargs abstol is preserved when call omits it" begin
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
