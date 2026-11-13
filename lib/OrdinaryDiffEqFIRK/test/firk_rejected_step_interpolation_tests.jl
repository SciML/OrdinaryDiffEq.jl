using OrdinaryDiffEqFIRK, OrdinaryDiffEqCore, Test

@testset "RadauIIA5 interpolation after error rejection" begin
    growth(u, p, t) = u
    growth!(du, u, p, t) = (du .= u)
    for f in (growth, growth!), extrapolant in (:dense, :constant)
        @testset "in-place=$(f === growth!), extrapolant=$extrapolant" begin
            integrator = init(
                ODEProblem(f, [1.0], (0.0, 2.0)), RadauIIA5(; extrapolant);
                dt = 1.0, reltol = 1.0e-12, abstol = 1.0e-12
            )
            integrator.iter = 1
            OrdinaryDiffEqCore.perform_step!(integrator, integrator.cache)
            @test !integrator.force_stepfail
            @test OrdinaryDiffEqCore.get_EEst(integrator) > 1
            endpoint = OrdinaryDiffEqCore.ode_interpolant(
                1.0, integrator.dt, integrator.uprev, integrator.u,
                integrator.k, integrator.cache, nothing, Val{0},
                integrator.differential_vars
            )
            @test endpoint == integrator.u
        end
    end
end
