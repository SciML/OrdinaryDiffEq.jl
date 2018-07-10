using OrdinaryDiffEq, DiffEqDevTools, Test, Random
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear

srand(100)
linear_bigαN = big"0.5"
f_linearbig = (u,p,t) -> (linear_bigαN*u)
f_2dlinearbig = (du,u,p,t) -> (du.=linear_bigαN*u)
(f::typeof(f_linearbig))(::Type{Val{:analytic}},u0,p,t) = u0*exp(linear_bigαN*t)
(f::typeof(f_2dlinearbig))(::Type{Val{:analytic}},u0,p,t) = u0*exp.(linear_bigαN*t)
probArr = [ODEProblem(f_linearbig, big"0.5", (0,1.)),
           ODEProblem(f_2dlinearbig, big.(rand(4,2)), (0,1.)),]
testTol = 0.2
dts = 1 .//(2 .^(10:-1:4))

@testset "Nordsieck Convergence Tests" begin
  for i in eachindex(probArr)
    sim = test_convergence(dts,probArr[i],AN5())
    @test abs(sim.𝒪est[:final]-5) < testTol
    @test abs(sim.𝒪est[:l2]-5) < testTol
    @test abs(sim.𝒪est[:l∞]-5) < testTol
  end
end

probArr = [prob_ode_linear,
           prob_ode_2Dlinear]
@testset "Nordsieck Adaptivity Tests: AN5" begin
  for i in eachindex(probArr)
    prob = probArr[i]
    sol = solve(prob, AN5(), reltol=1e-6)
    @test length(sol.t) < 11
    exact = prob.f(Val{:analytic}, prob.u0, prob.p, prob.tspan[end])
    @test Float64(norm(exact-sol[end])) < 1e-5
  end
end

@testset "Nordsieck Adaptivity Tests: JVODE" begin
  for i in eachindex(probArr)
    prob = probArr[i]
    sol = solve(prob, JVODE_Adams(), reltol=1e-4, abstol=1e-7)
    @test length(sol.t) < 22
    exact = prob.f(Val{:analytic}, prob.u0, prob.p, prob.tspan[end])
    @test norm(exact - sol[end], Inf) < 3e-3
  end
end
