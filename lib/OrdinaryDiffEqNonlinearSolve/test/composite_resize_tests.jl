using ADTypes: AutoFiniteDiff, AutoForwardDiff
using OrdinaryDiffEqBDF
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqCore: CompositeAlgorithm
using Test

function composite_resize_rhs!(du, u, p, t)
    fill!(du, 0)
    du[1] = sum(@view u[2:end])
    return nothing
end

@testset "Composite nonlinear resize algorithm context" begin
    algorithms = (
        (FBDF(), Tsit5()),
        (Tsit5(), FBDF()),
        (FBDF(autodiff = AutoForwardDiff()), FBDF(autodiff = AutoFiniteDiff())),
        (FBDF(autodiff = AutoFiniteDiff()), FBDF(autodiff = AutoForwardDiff())),
    )
    @testset "case=$case, active=$initial, $operation" for (case, components) in enumerate(algorithms),
            initial in 1:2, operation in (:resize, :delete, :add)
        active = Ref(initial)
        alg = CompositeAlgorithm(components, _ -> active[])
        integrator = init(
            ODEProblem(composite_resize_rhs!, [1.0, 2.0, 3.0, 4.0], (0.0, 1.0)),
            alg; adaptive = false, dt = 1 / 64
        )
        step!(integrator)
        before = copy(integrator.u)
        if operation == :delete
            deleteat!(integrator, [2])
            expected = before[[1, 3, 4]]
        else
            operation == :resize ? resize!(integrator, 5) : addat!(integrator, [2])
            integrator.u[5] = 5.0
            integrator.uprev[5] = 5.0
            expected = vcat(before, 5.0)
        end
        @test integrator.u == expected
        for cache in integrator.cache.caches
            if hasproperty(cache, :nlsolver)
                @test length(cache.nlsolver.z) == length(expected)
                @test length(cache.nlsolver.cache.du1) == length(expected)
            end
        end
        derivative_discontinuity!(integrator, true)
        for selection in (3 - initial, initial)
            active[] = selection
            step!(integrator)
            expected[1] += sum(@view expected[2:end]) / 64
            @test integrator.u ≈ expected rtol = 1.0e-10 atol = 1.0e-12
            @test integrator.cache.current == selection
        end
    end
end
