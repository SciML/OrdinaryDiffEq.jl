using OrdinaryDiffEqFIRK, DiffEqDevTools, Test, LinearAlgebra
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, van

testTol = 0.3

for prob in [prob_ode_linear, prob_ode_2Dlinear]
    sim21 = test_convergence(1 .// 2 .^ (6:-1:3), prob, RadauIIA5())
    @test sim21.𝒪est[:final]≈5 atol=testTol
end

sim21 = test_convergence(1 ./ 2 .^ (2.5:-1:0.5), prob_ode_linear, RadauIIA9())
@test sim21.𝒪est[:final]≈8.5 atol=testTol

sim21 = test_convergence(1 ./ 2 .^ (2.5:-1:0.5), prob_ode_2Dlinear, RadauIIA9())
@test sim21.𝒪est[:final]≈8.5 atol=testTol

prob_ode_linear_big = remake(prob_ode_linear, u0 = big.(prob_ode_linear.u0), tspan = big.(prob_ode_linear.tspan))
prob_ode_2Dlinear_big = remake(prob_ode_2Dlinear, u0 = big.(prob_ode_2Dlinear.u0), tspan = big.(prob_ode_2Dlinear.tspan))

for i in [3, 5, 7, 9], prob in [prob_ode_linear_big, prob_ode_2Dlinear_big]
    dts = 1 ./ 2 .^ (4.25:-1:0.25)
    sim21 = test_convergence(dts, prob, AdaptiveRadau(min_num_stages = i, max_num_stages = i))
    @test sim21.𝒪est[:final]≈ (2 * i - 1) atol=testTol
end

solve(prob_ode_linear, AdaptiveRadau(min_num_stages = 3, max_num_stages = 3))
solve(prob_ode_linear, RadauIIA5())

using OrdinaryDiffEq, StaticArrays, LinearSolve, ParameterizedFunctions

hires = @ode_def Hires begin
    dy1 = -1.71 * y1 + 0.43 * y2 + 8.32 * y3 + 0.0007
    dy2 = 1.71 * y1 - 8.75 * y2
    dy3 = -10.03 * y3 + 0.43 * y4 + 0.035 * y5
    dy4 = 8.32 * y2 + 1.71 * y3 - 1.12 * y4
    dy5 = -1.745 * y5 + 0.43 * y6 + 0.43 * y7
    dy6 = -280.0 * y6 * y8 + 0.69 * y4 + 1.71 * y5 - 0.43 * y6 + 0.69 * y7
    dy7 = 280.0 * y6 * y8 - 1.81 * y7
    dy8 = -280.0 * y6 * y8 + 1.81 * y7
end

u0 = SA[1, 0, 0, 0, 0, 0, 0, 0.0057]
probiip = ODEProblem{true}(hires, Vector(u0), (0.0, 10.0))
proboop = ODEProblem{false}(hires, Vector(u0), (0.0, 10.0))
proboop = ODEProblem{false}(hires, u0, (0.0, 10.0))

#=@btime =# sol = solve(proboop, AdaptiveRadau(min_num_stages = 3, max_num_stages = 7), reltol = 1e-1)
#=@btime =# sol = solve(proboop, RadauIIA5())

function rober!(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    du[3] = k₂ * y₂^2
    nothing
end
prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1e5), [0.04, 3e7, 1e4])
sol = solve(prob, AdaptiveRadau(min_num_stages = 3, max_num_stages =7))
sol2 = solve(prob, RadauIIA5())

using BenchmarkTools
#oop
@btime solve(prob_ode_linear, RadauIIA5(), adaptive = false, dt = 1e-2)
@btime solve(prob_ode_linear, AdaptiveRadau(num_stages = 3), adaptive = false, dt = 1e-2)

#ip
@btime solve(prob_ode_2Dlinear, RadauIIA5(), adaptive = false, dt = 1e-2)
@btime solve(prob_ode_2Dlinear, AdaptiveRadau(num_stages = 3), adaptive = false, dt = 1e-2)

VSCodeServer.@profview solve(prob_ode_linear, AdaptiveRadau(num_stages = 3), adaptive = false, dt = 1e-5)
VSCodeServer.@profview solve(prob_ode_linear, RadauIIA5(), adaptive = false, dt = 1e-5)


solve(prob_ode_linear, AdaptiveRadau(num_stages = 3), adaptive = false, dt = 1e-2)


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
