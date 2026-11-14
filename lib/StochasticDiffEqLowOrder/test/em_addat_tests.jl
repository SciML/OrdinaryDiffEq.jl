using StochasticDiffEqLowOrder
using SciMLBase: addat!
using Test

f!(du, u, p, t) = (du .= -0.5 .* u; nothing)
g!(du, u, p, t) = (du .= 0.1 .* u; nothing)

new_integrator(u0 = [1.0, 2.0]) = init(SDEProblem(f!, g!, u0, (0.0, 1.0)), EM(); dt = 1 / 64, seed = 1)

@testset "EM steps the added component from its value" begin
    integ = new_integrator()
    addat!(integ, 3:3, [1.0])
    @test integ.u[3] == integ.uprev[3] == 1.0
    step!(integ)
    # One EM step: u = uprev + dt * f(uprev) + g(uprev) * dW
    @test abs(integ.u[3] - (1.0 - 0.5 / 64)) < 0.1
end

@testset "EM keeps a zero added component at zero" begin
    integ = new_integrator([1.0, 2.0, 3.0])
    step!(integ)
    # Shrinking keeps a Vector's contents, so without a fill the old value comes back.
    deleteat!(integ, 3)
    addat!(integ, 3:3)
    @test integ.u[3] == integ.uprev[3] == 0.0
    solve!(integ)
    @test integ.u[3] == 0.0
end
