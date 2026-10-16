using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit
using LinearSolve: KrylovJL_GMRES
using SciMLBase: FullSpecialize, successful_retcode
using LinearAlgebra
using SparseArrays
using Test

#=
An implicit integrator using a Krylov linear solver must stay round-off
reproducible: the accepted step sequence may depend on `(W, b)`, but not on the
history of previous linear solves. A warm start that reuses an unpredictive
previous solution breaks that — the returned iterate is only accurate to
`reltol`, so it carries the whole history forward, and the step-size controller
then amplifies last-ulp input differences without bound. That is what
SciML/OrdinaryDiffEq.jl#4034 measured: perturbing `u0[1]` by one ulp moved the
`Hairer4` solution of this index-1 DAE by 3.7e-03.

The gate is deliberately loose relative to the values it protects. Different
BLAS kernels can perturb the nonlinear solve enough to select different step
sequences, and the drift is then set by the dense-output error of the algebraic
variables `u[3:4]` (1e-04 to 6e-04 here, against 1e-07 for `u[1:2]` and at the
step points), not by the linear solves. `Hairer4` drifts by 1.5e-04 this way
while its direct-factorization control also crosses 1e-05, and `Hairer42`, whose
smoothed error estimate settles into an accept-accept-reject cycle on this DAE,
drifts by up to 8.3e-05 with either GMRES or a direct factorization. The 5e-04
and 2.5e-04 bounds remain sevenfold below the 3.7e-03 `Hairer4` regression and
more than twofold below the 5.7e-04 `Hairer42` regression.
=#

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

@testset "Newton-Krylov round-off reproducibility ($(nameof(alg)))" for (alg, drift_limit) in (
        (Hairer4, 5.0e-4), (Hairer42, 2.5e-4),
    )
    a = solve(
        dae_prob(U0), alg(linsolve = KrylovJL_GMRES());
        reltol = 1.0e-8, abstol = 1.0e-6, maxiters = 10_000
    )
    b = solve(
        dae_prob(U0_PERTURBED), alg(linsolve = KrylovJL_GMRES());
        reltol = 1.0e-8, abstol = 1.0e-6, maxiters = 10_000
    )
    @test successful_retcode(a)
    @test successful_retcode(b)
    drift = maximum(
        maximum(abs.(a(t) .- b(t)))
            for t in DAE_TSPAN[1]:0.1:DAE_TSPAN[2]
    )
    @test drift < drift_limit
end
