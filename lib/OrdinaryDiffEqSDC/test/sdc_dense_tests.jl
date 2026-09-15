using OrdinaryDiffEqSDC
using SciMLBase
using Test

rotation!(du, u, p, t) = (du[1] = -u[2]; du[2] = u[1]; nothing)
const DENSE_ROTATION = ODEProblem(rotation!, [1.0, 0.0], (0.0, 4.0))
rotation_exact(t, ::Type{Val{0}}) = [cos(t), sin(t)]
rotation_exact(t, ::Type{Val{1}}) = [-sin(t), cos(t)]

function dense_order(M, dts, deriv)
    alg = SDC(num_nodes = M, num_sweeps = 2M)
    ts = range(DENSE_ROTATION.tspan..., length = 2001)
    errors = map(dts) do dt
        sol = solve(DENSE_ROTATION, alg; dt = dt, adaptive = false)
        sqrt(sum(sum(abs2, sol(t, deriv) .- rotation_exact(t, deriv)) for t in ts) / length(ts))
    end
    return log2(errors[1] / errors[2]) / log2(dts[1] / dts[2])
end

@testset "SDC dense output is the collocation polynomial" begin
    # With converged sweeps the polynomial through the node values has order M + 1
    # and its derivative order M, where cubic Hermite would cap both at 4 and 3.
    for (M, dts) in ((3, (0.2, 0.1)), (4, (0.4, 0.2)), (5, (0.4, 0.2)))
        @test dense_order(M, dts, Val{0}) ≈ M + 1 atol = 0.3
        @test dense_order(M, dts, Val{1}) ≈ M atol = 0.3
    end
end

@testset "SDC dense output interface" begin
    sol = solve(DENSE_ROTATION, SDC(num_nodes = 4, num_sweeps = 8); abstol = 1.0e-10, reltol = 1.0e-10)
    t = 1.2345
    out = zeros(2)
    sol(out, t)
    @test out == sol(t)
    @test sol(t; idxs = 2) == sol(t)[2]
    @test all(maximum(abs.(sol(prevfloat(sol.t[i])) .- sol.u[i])) < 1.0e-12 for i in 2:length(sol.t))
    @test SciMLBase.interp_summary(sol.interp) == "collocation polynomial"

    decay = solve(
        ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0)), SDC(num_nodes = 3, num_sweeps = 6);
        dt = 0.25, adaptive = false
    )
    @test abs(decay(0.3, Val{1}) + exp(-0.3)) < 1.0e-4

    # One Radau node makes the polynomial linear between the step values.
    single = solve(
        ODEProblem(rotation!, [1.0, 0.0], (0.0, 1.0)), SDC(num_nodes = 1, num_sweeps = 2);
        dt = 0.1, adaptive = false
    )
    i = searchsortedlast(single.t, 0.55)
    θ = (0.55 - single.t[i]) / (single.t[i + 1] - single.t[i])
    @test single(0.55) ≈ single.u[i] + θ * (single.u[i + 1] - single.u[i]) rtol = 1.0e-12
end
