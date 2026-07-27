using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit
using LinearSolve: KrylovJL_GMRES
using SciMLBase: FullSpecialize, successful_retcode
using LinearAlgebra
using SparseArrays
using Test

# A non-predictive warm start can carry Krylov error into adaptive step selection
# and amplify a one-ulp input perturbation. The threshold leaves ample margin for
# ordinary step-sequence sensitivity while rejecting macroscopic drift.

dae!(du, u, p, t) = mul!(du, p, u)

const DAE_P = [
    -1.0 0.0 0.0 0.0
    1.0 -0.5 0.0 0.0
    1.0 1.0 -1.0 0.0
    -1.0 1.0 0.0 -1.0
]
const DAE_MM = Diagonal([1.0, 1.0, 0.0, 0.0])
const DAE_JP = sparse(map(x -> iszero(x) ? 0.0 : 1.0, DAE_P))
const DAE_TSPAN = (0.0, 5.0)

dae_prob(u0) = ODEProblem(
    ODEFunction{true, FullSpecialize}(
        dae!; mass_matrix = DAE_MM, jac_prototype = DAE_JP
    ),
    u0, DAE_TSPAN, DAE_P; initializealg = BrownFullBasicInit()
)

const U0 = [1.0, 1.0, 0.5, 0.5]
const U0_PERTURBED = let u = copy(U0)
    u[1] = nextfloat(u[1])
    u
end

@testset "Newton-Krylov round-off reproducibility ($(nameof(alg)))" for alg in (
        Hairer4, Hairer42,
    )
    a = solve(dae_prob(U0), alg(linsolve = KrylovJL_GMRES()); maxiters = 10_000)
    b = solve(dae_prob(U0_PERTURBED), alg(linsolve = KrylovJL_GMRES()); maxiters = 10_000)
    @test successful_retcode(a)
    @test successful_retcode(b)
    drift = maximum(
        maximum(abs.(a(t) .- b(t)))
            for t in DAE_TSPAN[1]:0.1:DAE_TSPAN[2]
    )
    @test drift < 1.0e-5
end
