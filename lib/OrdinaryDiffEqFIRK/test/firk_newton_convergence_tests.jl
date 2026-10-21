using OrdinaryDiffEqFIRK, Test
using SciMLBase: successful_retcode

@testset "FIRK Newton convergence with a stale contraction estimate" begin
    function switched!(du, u, p, t)
        e = exp(-t)
        for i in eachindex(u)
            du[i] = -10i * (u[i] - e) - e
            if t > 0.2 && i < length(u)
                du[i] += 3 * (u[i + 1] - e)
            end
        end
    end
    function switched(u, p, t)
        du = similar(u)
        switched!(du, u, p, t)
        return du
    end

    stationary!(du, u, p, t) = fill!(du, 0)
    stationary(u, p, t) = zero(u)

    for alg in (
                RadauIIA3(), RadauIIA5(), RadauIIA9(),
                AdaptiveRadau(min_order = 9, max_order = 9),
                GaussLegendre(num_stages = 3),
            ), (f, steady) in ((switched!, stationary!), (switched, stationary))
        @testset "$(nameof(typeof(alg))) in-place=$(f === switched!)" begin
            prob = ODEProblem(f, ones(8), (0.0, 0.5))
            integrator = init(
                prob, alg; reltol = 1.0e-10, abstol = 1.0e-12,
                dt = 0.01, tstops = [0.2]
            )
            while integrator.t < 0.2
                step!(integrator)
            end
            integrator.cache.ηold = eps(Float64)
            step!(integrator)
            @test integrator.cache.iter >= 2
            solve!(integrator)
            @test successful_retcode(integrator.sol)
            @test maximum(abs.(integrator.u .- exp(-0.5))) < 1.0e-8

            equilibrium = init(ODEProblem(steady, ones(8), (0.0, 0.5)), alg; dt = 0.01)
            step!(equilibrium)
            @test equilibrium.cache.iter == 1
            @test equilibrium.u == ones(8)
        end
    end
end
