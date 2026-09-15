using OrdinaryDiffEqSDC
using SciMLBase
using Test

rotation!(du, u, p, t) = (du[1] = -u[2]; du[2] = u[1]; nothing)
const DENSE_ROTATION = ODEProblem(rotation!, [1.0, 0.0], (0.0, 4.0))
rotation_exact(t, ::Type{Val{0}}) = [cos(t), sin(t)]
rotation_exact(t, ::Type{Val{1}}) = [-sin(t), cos(t)]

function dense_orders(M, dts)
    alg = SDC(num_nodes = M, num_sweeps = 2M)
    ts = range(DENSE_ROTATION.tspan..., length = 2001)
    errors = map(dts) do dt
        sol = solve(DENSE_ROTATION, alg; dt = dt, adaptive = false)
        return map((Val{0}, Val{1})) do deriv
            sqrt(
                sum(sum(abs2, sol(t, deriv) .- rotation_exact(t, deriv)) for t in ts) /
                    length(ts)
            )
        end
    end
    return map(
        (first, last)
    ) do part
        log2(part(errors[1]) / part(errors[2])) / log2(dts[1] / dts[2])
    end
end

@testset "SDC dense output is the collocation polynomial" begin
    # With converged sweeps the polynomial through the node values has order M + 1
    # and its derivative order M, where cubic Hermite would cap both at 4 and 3.
    for (M, dts) in ((3, (0.2, 0.1)), (4, (0.4, 0.2)), (5, (0.4, 0.2)))
        value_order, derivative_order = dense_orders(M, dts)
        @test value_order ≈ M + 1 atol = 0.3
        @test derivative_order ≈ M atol = 0.3
    end
end

@testset "SDC dense output falls back when the node rates do not describe the step" begin
    decay = ODEProblem((du, u, p, t) -> (du[1] = -u[1]; nothing), [1.0], (0.0, 1.0))
    alg = SDC(num_nodes = 3, num_sweeps = 6)

    # A ContinuousCallback shortens the step it fires in, which leaves the node rates
    # describing an interval the solution no longer covers.
    cb = ContinuousCallback((u, t, i) -> u[1] - 0.5, i -> nothing)
    sol = solve(decay, alg; dt = 0.5, adaptive = false, callback = cb, saveat = 0.0:0.05:1.0)
    @test maximum(abs(sol.u[i][1] - exp(-sol.t[i])) for i in eachindex(sol.t)) < 1.0e-3

    event = solve(decay, alg; dt = 0.5, adaptive = false, callback = cb)
    i = findfirst(u -> isapprox(u[1], 0.5; atol = 1.0e-8), event.u)
    @test event(prevfloat(event.t[i]))[1] ≈ event.u[i][1] rtol = 1.0e-8

    # `LastNode` ends the step on the last node rather than on the polynomial.
    last_node = solve(
        ODEProblem(rotation!, [1.0, 0.0], (0.0, 4.0)),
        SDC(num_nodes = 3, num_sweeps = 1, step_update = SDCStepUpdate.LastNode);
        dt = 0.2, adaptive = false
    )
    @test maximum(
        maximum(abs.(last_node(prevfloat(last_node.t[i])) .- last_node.u[i]))
            for i in 2:length(last_node.t)
    ) < 1.0e-10
end

@testset "SDC hands the interpolant derivatives" begin
    prob = ODEProblem(rotation!, [1.0, 0.0], (0.0, 4.0))
    integrator = init(prob, SDC(num_nodes = 3, num_sweeps = 3); dt = 0.2, adaptive = false)
    @test SciMLBase.get_du(integrator) ≈ [0.0, 1.0] atol = 1.0e-12

    sol = solve(prob, SDC(num_nodes = 3, num_sweeps = 6); dt = 0.2, adaptive = false)
    # `k[1]` and `k[2]` are the derivatives at the ends of the step, in rate units.
    @test sol.k[2][1] ≈ [-sol.u[1][2], sol.u[1][1]] rtol = 1.0e-3
    @test sol.k[2][2] ≈ [-sol.u[2][2], sol.u[2][1]] rtol = 1.0e-3
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
