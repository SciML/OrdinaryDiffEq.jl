using OrdinaryDiffEqStabilizedRK, Test
import OrdinaryDiffEqCore
using SciMLBase: DEIntegrator, ReturnCode

const STABILIZED_ALGS = (
    ROCK2, ROCK4, RKC, RKMC2, ESERK4, ESERK5, SERK2, TSRKC2, TSRKC3,
    RKL1, RKL2, RKG1, RKG2,
)

@testset "has_stage_limiter opt-in" begin
    for A in STABILIZED_ALGS
        @test OrdinaryDiffEqCore.has_stage_limiter(A())
    end
end

# A clamp applied by the stage limiter must be visible from inside `f` for every
# state handed to it during a step, and the accepted `u` must satisfy it too.
@testset "every stage value is limited before f sees it" begin
    eigen_est = integrator -> integrator.eigen_est = 500.0
    floor = 0.9
    for A in STABILIZED_ALGS, adaptive in (false, true)
        alg = A(; eigen_est)
        violations = Ref(0)
        calls = Ref(0)
        got_integrator = Ref(true)
        got_uprev = Ref(false)
        function f!(du, u, p, t)
            any(<(floor), u) && (violations[] += 1)
            @. du = -10 * u
            return nothing
        end
        function limiter!(u, integrator, p, t)
            calls[] += 1
            integrator isa DEIntegrator || (got_integrator[] = false)
            u === integrator.uprev && (got_uprev[] = true)
            @. u = max(u, floor)
            return nothing
        end
        prob = ODEProblem(f!, [1.0, 1.0], (0.0, 1.0))
        sol = solve(prob, alg; dt = 0.1, adaptive, stage_limiter = limiter!)
        @test sol.retcode == ReturnCode.Success
        @test calls[] > 0
        @test got_integrator[]
        @test !got_uprev[]
        @test violations[] == 0
        @test all(u -> all(>=(floor), u), sol.u)

        unlimited = solve(prob, alg; dt = 0.1, adaptive)
        @test unlimited.u[end][1] < floor
    end
end
