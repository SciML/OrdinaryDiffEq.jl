using OrdinaryDiffEqFIRK, DiffEqDevTools, Test, LinearAlgebra
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, prob_ode_vanderpol

testTol = 0.5

for prob in [prob_ode_linear, prob_ode_2Dlinear]
    sim21 = test_convergence(1 .// 2 .^ (6:-1:3), prob, RadauIIA5(), dense_errors = true)
    @test sim21.𝒪est[:final]≈5 atol=testTol
    @test sim21.𝒪est[:L2]≈4 atol=testTol
end

sim21 = test_convergence(
    1 ./ 2 .^ (2.5:-1:0.5), prob_ode_linear, RadauIIA9(), dense_errors = true)
@test sim21.𝒪est[:final]≈8 atol=testTol
@test sim21.𝒪est[:L2]≈6 atol=testTol

sim21 = test_convergence(
    1 ./ 2 .^ (2.5:-1:0.5), prob_ode_2Dlinear, RadauIIA9(), dense_errors = true)
@test sim21.𝒪est[:final]≈8 atol=testTol
@test sim21.𝒪est[:L2]≈6 atol=testTol

using GenericSchur

prob_ode_linear_big = remake(
    prob_ode_linear, u0 = big.(prob_ode_linear.u0), tspan = big.(prob_ode_linear.tspan))
prob_ode_2Dlinear_big = remake(prob_ode_2Dlinear, u0 = big.(prob_ode_2Dlinear.u0),
    tspan = big.(prob_ode_2Dlinear.tspan))

#non-threaded tests
for i in [5, 9, 13, 17, 21, 25], prob in [prob_ode_linear_big, prob_ode_2Dlinear_big]
    dts = 1 ./ 2 .^ (4.25:-1:0.25)
    local sim21 = test_convergence(
        dts, prob, AdaptiveRadau(min_order = i, max_order = i), dense_errors = true)
    @test sim21.𝒪est[:final]≈i atol=testTol
    @test sim21.𝒪est[:L2]≈((i + 3) ÷ 2) atol=testTol
end

#threaded tests
using OrdinaryDiffEqCore
for i in [5, 9, 13, 17, 21, 25], prob in [prob_ode_linear_big, prob_ode_2Dlinear_big]
    dts = 1 ./ 2 .^ (4.25:-1:0.25)
    local sim21 = test_convergence(dts,
        prob,
        AdaptiveRadau(min_order = i, max_order = i,
            threading = OrdinaryDiffEqCore.PolyesterThreads()))
    @test sim21.𝒪est[:final]≈i atol=testTol
end

sys = prob_ode_vanderpol.f.sys

# test adaptivity
for iip in (true, false)
    vanstiff = ODEProblem{iip}(
        sys, [sys.y => 0, sys.x => sqrt(3), sys.μ => 1e6], (0.0, 1.0))
    sol = solve(vanstiff, RadauIIA5())
    if iip
        @test sol.stats.naccept + sol.stats.nreject > sol.stats.njacs # J reuse
        @test sol.stats.njacs < sol.stats.nw # W reuse
    end
    @test length(sol) < 150
    @test SciMLBase.successful_retcode(sol)
    sol_temp = solve(remake(vanstiff, p = [sys.μ => 1e7]), RadauIIA5())
    @test length(sol_temp) < 150
    @test SciMLBase.successful_retcode(sol_temp)
    sol_temp2 = solve(
        remake(vanstiff, p = [sys.μ => 1e7]), reltol = [1e-6, 1e-4], RadauIIA5())
    @test length(sol_temp2) < 180
    @test SciMLBase.successful_retcode(sol_temp2)
    sol_temp3 = solve(remake(vanstiff, p = [sys.μ => 1e7]), RadauIIA5(), reltol = 1e-9,
        abstol = 1e-9)
    @test length(sol_temp3) < 970
    @test SciMLBase.successful_retcode(sol_temp3)
    sol_temp4 = solve(remake(vanstiff, p = [sys.μ => 1e9]), RadauIIA5())
    @test length(sol_temp4) < 170
    @test SciMLBase.successful_retcode(sol_temp4)
    sol_temp5 = solve(remake(vanstiff, p = [sys.μ => 1e10]), RadauIIA5())
    @test length(sol_temp5) < 190
    @test SciMLBase.successful_retcode(sol_temp5)
end

##Tests for RadauIIA3
for prob in [prob_ode_linear, prob_ode_2Dlinear]
    dts = 1 ./ 2 .^ (8:-1:1)
    sim = test_convergence(dts, prob, RadauIIA3(), dense_errors = true)
    @test sim.𝒪est[:final]≈3 atol=0.25
    @test sim.𝒪est[:L2]≈3 atol=0.25
end

# test adaptivity
for iip in (true, false)
    vanstiff = ODEProblem{iip}(
        sys, [sys.y => 0, sys.x => sqrt(3), sys.μ => 1e6], (0.0, 1.0))
    sol = solve(vanstiff, RadauIIA3())
    if iip
        @test sol.stats.naccept + sol.stats.nreject > sol.stats.njacs # J reuse
        @test sol.stats.njacs < sol.stats.nw # W reuse
    end
    @test length(sol) < 5000 # the error estimate is not very good
    @test SciMLBase.successful_retcode(sol)
end
