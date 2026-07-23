using OrdinaryDiffEqDifferentiation: default_krylov_warm_start
using LinearSolve
using Test

@testset "default_krylov_warm_start resolves Auto for the Newton path" begin
    # Non-Krylov solvers are always returned unchanged.
    @test default_krylov_warm_start(LUFactorization()) isa LUFactorization

    if isdefined(LinearSolve, :WarmStart)
        # Auto (the KrylovJL default) resolves to Hegedus on the Newton path.
        @test default_krylov_warm_start(KrylovJL_GMRES()).warm_start ===
            LinearSolve.WarmStart.Hegedus
        # Explicit choices are respected, never overridden.
        for ws in (
                LinearSolve.WarmStart.None, LinearSolve.WarmStart.Previous,
                LinearSolve.WarmStart.Hegedus,
            )
            @test default_krylov_warm_start(KrylovJL_GMRES(warm_start = ws)).warm_start === ws
        end
        # Other fields survive the remake.
        precs = (A, p) -> (I, I)
        resolved = default_krylov_warm_start(KrylovJL_GMRES(; precs))
        @test resolved.warm_start === LinearSolve.WarmStart.Hegedus
        @test resolved.precs === precs
    else
        # LinearSolve too old to define WarmStart: strictly a no-op.
        k = KrylovJL_GMRES()
        @test default_krylov_warm_start(k) === k
    end
end
