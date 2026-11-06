using OrdinaryDiffEq, Test, ADTypes
using OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK
import ODEProblemLibrary: prob_ode_vanderpol
using ForwardDiff: Dual, Tag
using LinearAlgebra: mul!
import OrdinaryDiffEqCore

# Create Van der Pol problem with same structure as the new ODEProblemLibrary implementation
# New implementation uses: u[1] = x, u[2] = y, p[1] = μ
# Van der Pol equations: dx/dt = y, dy/dt = μ * ((1 - x^2) * y - x)
# Initial conditions: [x, y] = [2.0, 0] (matching the original [sys.x => 2.0, sys.y => 0])
function __van(du, u, p, t)
    x, y = u[1], u[2]
    μ = p[1]
    du[1] = y                           # dx/dt = y
    du[2] = μ * ((1 - x^2) * y - x)     # dy/dt = μ * ((1 - x^2) * y - x)
    return
end
prob1 = ODEProblem(__van, [2.0, 0.0], (0.0, 6), [inv(0.003)])
prob2 = ODEProblem(__van, [2.0, 0.0], (0.0, 6), [inv(0.003)])
# out-of-place test
function _van(u, p, t)
    x, y = u[1], u[2]
    μ = p[1]
    return [
        y,                           # dx/dt = y
        μ * ((1 - x^2) * y - x),
    ]     # dy/dt = μ * ((1 - x^2) * y - x)
end
prob3 = ODEProblem(_van, [2.0, 0.0], (0.0, 6), [inv(0.003)])
probArr = [prob1, prob2, prob3]

for prob in [prob2, prob3], u0 in [prob.u0, Dual.(prob.u0, prob.u0)]
    prob′ = remake(prob3; u0)
    @test_nowarn solve(prob′, AutoTsit5(Rosenbrock23(autodiff = AutoFiniteDiff())))
end

# Test if switching back and forth
is_switching_fb(sol) = all(i -> count(isequal(i), sol.alg_choice[2:end]) > 5, (1, 2))
for (i, prob) in enumerate(probArr)
    println(i)
    sol = solve(
        prob, AutoTsit5(Rosenbrock23(autodiff = AutoFiniteDiff())),
        maxiters = 1000
    )
    @test is_switching_fb(sol)
    alg = AutoTsit5(Rodas5(); maxstiffstep = 5, maxnonstiffstep = 5, stiffalgfirst = true)
    sol = solve(prob, alg, maxiters = 1000)
    sol2 = solve(prob, alg, maxiters = 1000)
    @test sol.t == sol2.t # test reinitialization
    @test length(sol.t) < 280
    @test SciMLBase.successful_retcode(sol)
    @test alg.algs[sol.alg_choice[1]] isa Rodas5
    i == 1 || @test is_switching_fb(sol) # fails due to eigenvalue estimate of J
    sol = solve(
        prob,
        AutoDP5(
            Rodas5(); maxstiffstep = 2, maxnonstiffstep = 2,
            stifftol = 11 // 10, nonstifftol = 9 / 10
        ),
        reltol = 1.0e-5, abstol = 1.0e-5, maxiters = 1000
    )
    @test length(sol.t) < 625
    @test SciMLBase.successful_retcode(sol)
    @test is_switching_fb(sol)

    sol = solve(
        prob, AutoVern6(Kvaerno3(); maxstiffstep = 4, maxnonstiffstep = 2),
        maxiters = 1000
    )
    @test length(sol.t) < 700
    @test SciMLBase.successful_retcode(sol)
    @test is_switching_fb(sol)
    sol = solve(
        prob, AutoVern7(Hairer42(); maxstiffstep = 4, maxnonstiffstep = 2),
        maxiters = 1000
    )
    @test length(sol.t) < 610
    @test_skip SciMLBase.successful_retcode(sol)
    @test is_switching_fb(sol)
    sol = solve(
        prob, AutoVern8(Rosenbrock23(); maxstiffstep = 4, maxnonstiffstep = 4),
        maxiters = 1000
    )
    @test length(sol.t) < 910
    @test SciMLBase.successful_retcode(sol)
    @test is_switching_fb(sol)
    sol = solve(
        prob, AutoVern9(KenCarp3(); maxstiffstep = 4, maxnonstiffstep = 1),
        maxiters = 1000
    )
    @test length(sol.t) < 570
    @test SciMLBase.successful_retcode(sol)
    @test is_switching_fb(sol)
    sol = solve(
        prob,
        AutoVern9(
            KenCarp3(autodiff = AutoFiniteDiff()); maxstiffstep = 4,
            maxnonstiffstep = 1
        ), maxiters = 1000
    )
    @test length(sol.t) < 570
    @test SciMLBase.successful_retcode(sol)
    @test is_switching_fb(sol)
end

# Stiffness detection must not depend on the partials of dual numbers (#4776):
# the forward sensitivity system has the same eigenvalues as the primal problem
@testset "Stiffness detection ignores dual partials" begin
    stiffalg = Rosenbrock23(autodiff = AutoFiniteDiff())
    algs = (
        AutoTsit5(stiffalg), AutoDP5(stiffalg), AutoVern6(stiffalg),
        AutoVern7(stiffalg), AutoVern8(stiffalg), AutoVern9(stiffalg),
    )
    seed = 1.0e3
    k1 = Dual{typeof(Tag(identity, Float64))}(0.25, seed)
    k2 = Dual{typeof(Tag(abs, typeof(k1)))}(k1, Dual{typeof(Tag(identity, Float64))}(seed, seed))
    # `|p| = 0.25` for real and complex `p`; `[one(p), zero(p)]` also covers #4878
    for k in (k1, k2), p in (k, complex(0.6, 0.8) * k), u0 in ([one(p)], [one(p), zero(p)])
        prob_iip = ODEProblem((du, u, p, t) -> (du .= -p .* u), u0, (0.0, 40.0), p)
        prob_oop = ODEProblem((u, p, t) -> -p .* u, u0, (0.0, 40.0), p)
        for prob in (prob_iip, prob_oop)
            for alg in algs
                sol = solve(prob, alg; abstol = 1.0e-12, reltol = 1.0e-8)
                @test SciMLBase.successful_retcode(sol)
                @test sol.stats.maxeig ≈ 0.25 rtol = 1.0e-2
                @test all(==(1), sol.alg_choice)
            end

            # The default solver switches via its own logic and does not record `maxeig`
            sol = solve(prob; abstol = 1.0e-12, reltol = 1.0e-8)
            @test SciMLBase.successful_retcode(sol)
            @test allequal(sol.alg_choice)
        end
    end
end

# Components with identical stage values must not give `0/0 = NaN`, which is treated as stiff (#4878)
@testset "Stiffness detection skips components with identical stage values" begin
    @test OrdinaryDiffEqCore.eigen_est_from_stages([1.0, 0.0], zeros(2), [4.0, 0.0], zeros(2)) == 0.25
    @test OrdinaryDiffEqCore.eigen_est_from_stages([1.0, 1.0], zeros(2), [4.0, 0.0], zeros(2)) == Inf
    x = Dual{typeof(Tag(identity, Float64))}(1.0, 1.0)
    z = [zero(complex(x))]
    est = OrdinaryDiffEqCore.eigen_est_from_stages([complex(0.6x, 0.8x)], z, [complex(4x)], z)
    @test est isa Float64
    @test est ≈ 0.25

    # ForwardDiff does not support complex in-place problems
    stiffalg = Rosenbrock23(autodiff = AutoFiniteDiff())
    algs = (
        AutoTsit5(stiffalg), AutoDP5(stiffalg), AutoVern6(stiffalg),
        AutoVern7(stiffalg), AutoVern8(stiffalg), AutoVern9(stiffalg),
    )
    # `u' = A u` with dominant eigenvalue `λ`
    linear_problems(λ) = (
        ([λ 0; 0 λ], [1, 0]), # identically zero component
        ([λ 0; 0 0], [1, 5]), # component in steady state
        ([-0.1 0; 0.1 λ], [0, 1]), # empty depot compartment
    )
    for λ in (-0.25, complex(-0.15, 0.2)), (A, u0) in linear_problems(λ)
        u0 = convert(Vector{eltype(A)}, u0)
        prob_iip = ODEProblem((du, u, A, t) -> mul!(du, A, u), u0, (0.0, 40.0), A)
        prob_oop = ODEProblem((u, A, t) -> A * u, u0, (0.0, 40.0), A)
        for prob in (prob_iip, prob_oop)
            for alg in algs
                sol = solve(prob, alg; abstol = 1.0e-12, reltol = 1.0e-8)
                @test SciMLBase.successful_retcode(sol)
                @test sol.stats.maxeig ≈ 0.25 rtol = 1.0e-2
                @test all(==(1), sol.alg_choice)
            end

            sol = solve(prob; abstol = 1.0e-12, reltol = 1.0e-8)
            @test SciMLBase.successful_retcode(sol)
            @test allequal(sol.alg_choice)
        end
    end
end
