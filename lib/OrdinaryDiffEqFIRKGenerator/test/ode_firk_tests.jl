using OrdinaryDiffEqFIRK, OrdinaryDiffEqFIRKGenerator, DiffEqDevTools, Test, LinearAlgebra
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear

testTol = 0.5

prob_ode_linear_big = remake(prob_ode_linear, u0 = big.(prob_ode_linear.u0), tspan = big.(prob_ode_linear.tspan))
prob_ode_2Dlinear_big = remake(prob_ode_2Dlinear, u0 = big.(prob_ode_2Dlinear.u0), tspan = big.(prob_ode_2Dlinear.tspan))

for i in [17, 21], prob in [prob_ode_linear_big, prob_ode_2Dlinear_big]
    dts = 1 ./ 2 .^ (4.25:-1:0.25)
    sim21 = test_convergence(dts, prob, AdaptiveRadau(min_order = i, max_order = i))
    @test sim21.𝒪est[:final]≈ i atol=testTol
end
