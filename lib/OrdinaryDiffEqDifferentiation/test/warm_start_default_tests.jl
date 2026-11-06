using OrdinaryDiffEqDifferentiation: default_krylov_warm_start
using LinearSolve
using Test

@testset "default_krylov_warm_start resolves Auto for the Newton path" begin
    # Non-Krylov solvers are always returned unchanged.
    @test default_krylov_warm_start(LUFactorization()) isa LUFactorization

    # Auto (the KrylovJL default) resolves to Hegedus on the Newton path. This is
    # sound only because the [compat] floor pins a LinearSolve that discards an
    # unpredictive Hegedus guess; without it warm starting amplifies round-off (#4034).
    expected = LinearSolve.WarmStart.Hegedus
    @test default_krylov_warm_start(KrylovJL_GMRES()).warm_start === expected

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
    @test resolved.warm_start === expected
    @test resolved.precs === precs
end
