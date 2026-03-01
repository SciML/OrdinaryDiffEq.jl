using OrdinaryDiffEq, Test

# Test that solution arrays release excess allocated memory after solve
# (https://github.com/SciML/DiffEqBase.jl/issues/1289)
@testset "Solution memory release after solve" begin
    function lorenz!(du, u, p, t)
        du[1] = 10.0 * (u[2] - u[1])
        du[2] = u[1] * (28.0 - u[3]) - u[2]
        du[3] = u[1] * u[2] - (8 / 3) * u[3]
    end

    u0 = [1.0, 0.0, 0.0]
    tspan = (0.0, 100.0)
    prob = ODEProblem(lorenz!, u0, tspan)
    sol = solve(prob, Tsit5())

    n = length(sol.t)
    # After solve, allocated memory should be close to actual data size.
    # sizehint! in postamble! releases excess capacity from array growth during integration.
    t_data_size = n * sizeof(eltype(sol.t)) + sizeof(sol.t)
    @test Base.summarysize(sol.t) < 1.5 * t_data_size
end
