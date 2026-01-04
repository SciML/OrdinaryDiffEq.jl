using OrdinaryDiffEq, Test, ADTypes
import ODEProblemLibrary: prob_ode_vanderpol
using ForwardDiff: Dual

# Create Van der Pol problem with same structure as the new ODEProblemLibrary implementation
# New implementation uses: u[1] = x, u[2] = y, p[1] = μ
# Van der Pol equations: dx/dt = y, dy/dt = μ * ((1 - x^2) * y - x)
# Initial conditions: [x, y] = [2.0, 0] (matching the original [sys.x => 2.0, sys.y => 0])
function __van(du, u, p, t)
    x, y = u[1], u[2]
    μ = p[1]
    du[1] = y                           # dx/dt = y
    return du[2] = μ * ((1 - x^2) * y - x)     # dy/dt = μ * ((1 - x^2) * y - x)
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
    prob′ = remake(prob3, u0 = u0)
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
