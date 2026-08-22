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

Both solves ask for a tolerance rather than taking the default. The linear solve
inherits `reltol`, and at the default 1e-3 a Krylov solve returns a Newton
direction accurate only to 1e-3, which cannot drive the stage to 1e-3: `dt`
collapses and both methods return `Unstable` before this measures anything
(SciML/OrdinaryDiffEq.jl#4293).

The gate is deliberately loose relative to the values it protects: at 1e-12 a
round-off-reproducible integration lands at 6.0e-07 (`Hairer4`) and 4.6e-07
(`Hairer42`), the step sequence itself being mildly chaotic, while the regression
it catches sits at 5.7e-04 (`Hairer42`) and 3.7e-03 (`Hairer4`).
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

@testset "Newton-Krylov round-off reproducibility ($(nameof(alg)))" for alg in (
        Hairer4, Hairer42,
    )
    a = solve(
        dae_prob(U0), alg(linsolve = KrylovJL_GMRES());
        reltol = 1.0e-12, abstol = 1.0e-12, maxiters = 10_000
    )
    b = solve(
        dae_prob(U0_PERTURBED), alg(linsolve = KrylovJL_GMRES());
        reltol = 1.0e-12, abstol = 1.0e-12, maxiters = 10_000
    )
    @test successful_retcode(a)
    @test successful_retcode(b)
    drift = maximum(
        maximum(abs.(a(t) .- b(t)))
            for t in DAE_TSPAN[1]:0.1:DAE_TSPAN[2]
    )
    @test drift < 1.0e-5
end
