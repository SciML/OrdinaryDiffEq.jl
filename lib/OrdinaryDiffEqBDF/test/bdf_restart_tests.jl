using OrdinaryDiffEqBDF, SciMLBase, Test

@testset "BDF history after a discontinuity" begin
    for (name, alg) in ("QNDF" => QNDF(), "QNDF2" => QNDF2(), "QBDF" => QBDF(), "QBDF2" => QBDF2()),
            inplace in (false, true),
            direction in (1.0, -1.0), adaptive in (false, true), modification in (:state, :rewind)
        @testset "$name inplace=$inplace direction=$direction adaptive=$adaptive modification=$modification" begin
            f = inplace ? ((du, u, p, t) -> (du .= p[1] .* u .+ p[2])) :
                ((u, p, t) -> p[1] .* u .+ p[2])
            prob = ODEProblem(f, [1.0], (0.0, 3direction), [1.0, 0.0])
            integrator = init(prob, alg; dt = direction / 8, dtmax = 1 / 8, adaptive, save_everystep = false, abstol = 1.0e-10, reltol = 1.0e-10)
            while direction * integrator.t < 1
                step!(integrator)
            end
            if modification == :rewind
                change_t_via_interpolation!(integrator, (integrator.tprev + integrator.t) / 2)
            end
            event_time = integrator.t
            slope = startswith(name, "QBDF") ? 1.0 : 0.0
            integrator.u .= 4
            integrator.p .= (0, slope)
            derivative_discontinuity!(integrator, true)
            for _ in 1:8
                step!(integrator)
                @test integrator.u ≈ [4 + slope * (integrator.t - event_time)]
            end
        end
    end
end
