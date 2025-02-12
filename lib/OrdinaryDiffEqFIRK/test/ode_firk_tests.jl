using OrdinaryDiffEqFIRK, DiffEqDevTools, Test, LinearAlgebra
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, van

testTol = 0.5

for prob in [prob_ode_linear, prob_ode_2Dlinear]
    sim21 = test_convergence(1 .// 2 .^ (6:-1:3), prob, RadauIIA5())
    @test sim21.𝒪est[:final]≈5 atol=testTol
end

sim21 = test_convergence(1 ./ 2 .^ (2.5:-1:0.5), prob_ode_linear, RadauIIA9())
@test sim21.𝒪est[:final]≈8 atol=testTol

sim21 = test_convergence(1 ./ 2 .^ (2.5:-1:0.5), prob_ode_2Dlinear, RadauIIA9())
@test sim21.𝒪est[:final]≈8 atol=testTol

using GenericSchur

prob_ode_linear_big = remake(
    prob_ode_linear, u0 = big.(prob_ode_linear.u0), tspan = big.(prob_ode_linear.tspan))
prob_ode_2Dlinear_big = remake(prob_ode_2Dlinear, u0 = big.(prob_ode_2Dlinear.u0),
    tspan = big.(prob_ode_2Dlinear.tspan))

#non-threaded tests
for i in [5, 9, 13, 17, 21, 25], prob in [prob_ode_linear_big, prob_ode_2Dlinear_big]
    dts = 1 ./ 2 .^ (4.25:-1:0.25)
    local sim21 = test_convergence(dts, prob, AdaptiveRadau(min_order = i, max_order = i))
    @test sim21.𝒪est[:final]≈i atol=testTol
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

# test adaptivity
for iip in (true, false)
    if iip
        vanstiff = ODEProblem{iip}(van, [0; sqrt(3)], (0.0, 1.0), 1e6)
    else
        vanstiff = ODEProblem{false}((u, p, t) -> van(u, p, t), [0; sqrt(3)], (0.0, 1.0),
            1e6)
    end
    sol = solve(vanstiff, RadauIIA5())
    if iip
        @test sol.stats.naccept + sol.stats.nreject > sol.stats.njacs # J reuse
        @test sol.stats.njacs < sol.stats.nw # W reuse
    end
    @test length(sol) < 150
    @test length(solve(remake(vanstiff, p = 1e7), RadauIIA5())) < 150
    @test length(solve(remake(vanstiff, p = 1e7), reltol = [1e-4, 1e-6], RadauIIA5())) < 170
    @test length(solve(remake(vanstiff, p = 1e7), RadauIIA5(), reltol = 1e-9,
        abstol = 1e-9)) < 870
    @test length(solve(remake(vanstiff, p = 1e9), RadauIIA5())) < 170
    @test length(solve(remake(vanstiff, p = 1e10), RadauIIA5())) < 190
end

##Tests for RadauIIA3
for prob in [prob_ode_linear, prob_ode_2Dlinear]
    dts = 1 ./ 2 .^ (8:-1:1)
    sim = test_convergence(dts, prob, RadauIIA3())
    @test sim.𝒪est[:final]≈3 atol=0.25
end

# test adaptivity
for iip in (true, false)
    if iip
        vanstiff = ODEProblem{iip}(van, [0; sqrt(3)], (0.0, 1.0), 1e6)
    else
        vanstiff = ODEProblem{false}((u, p, t) -> van(u, p, t), [0; sqrt(3)], (0.0, 1.0),
            1e6)
    end
    sol = solve(vanstiff, RadauIIA3())
    if iip
        @test sol.stats.naccept + sol.stats.nreject > sol.stats.njacs # J reuse
        @test sol.stats.njacs < sol.stats.nw # W reuse
    end
    @test length(sol) < 5000 # the error estimate is not very good
end
