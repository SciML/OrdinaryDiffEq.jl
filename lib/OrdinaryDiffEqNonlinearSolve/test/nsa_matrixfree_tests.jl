using OrdinaryDiffEqBDF, OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using NonlinearSolve, LinearSolve, LinearAlgebra, SparseArrays, ADTypes
using SparseConnectivityTracer, SciMLBase
using SciMLOperators: AbstractSciMLOperator
using Test

# With a Krylov inner linear solver the ODE builds a matrix-free `WOperator`; NonlinearSolveAlg
# now hands that operator to NonlinearSolve, which applies it via `mul!`. Previously the inner
# solve rebuilt a residual-derived AD Jacobian-vector product, which crashed under ForwardDiff
# (`No matching function wrapper`) because it differentiated the FunctionWrapped ODE residual.

const N = 8
const xyd = range(0, stop = 1, length = N)
lim(a, n) = a == n + 1 ? 1 : a == 0 ? n : a
function bruss!(du, u, p, t)
    A, B, al, dx = p
    al = al / dx^2
    return @inbounds for I in CartesianIndices((N, N))
        i, j = Tuple(I)
        ip, im, jp, jm = lim(i + 1, N), lim(i - 1, N), lim(j + 1, N), lim(j - 1, N)
        du[i, j, 1] = al * (
            u[im, j, 1] + u[ip, j, 1] + u[i, jp, 1] + u[i, jm, 1] -
                4u[i, j, 1]
        ) + B + u[i, j, 1]^2 * u[i, j, 2] - (A + 1) * u[i, j, 1]
        du[i, j, 2] = al * (
            u[im, j, 2] + u[ip, j, 2] + u[i, jp, 2] + u[i, jm, 2] -
                4u[i, j, 2]
        ) + A * u[i, j, 1] - u[i, j, 1]^2 * u[i, j, 2]
    end
end
u0 = ones(N, N, 2)
p = (3.4, 1.0, 10.0, step(xyd))
jsp = ADTypes.jacobian_sparsity(
    (du, u) -> bruss!(du, u, p, 0.0), similar(u0), u0,
    TracerSparsityDetector()
)
prob = ODEProblem(ODEFunction(bruss!; jac_prototype = float.(jsp)), u0, (0.0, 1.0), p)
refsol = solve(prob, Rodas5P(); abstol = 1.0e-10, reltol = 1.0e-10)

@testset "NSA + Krylov linsolve works for every inner autodiff (was a ForwardDiff crash)" begin
    for ad in (AutoFiniteDiff(), AutoForwardDiff())
        nsa = NonlinearSolveAlg(NewtonRaphson(; autodiff = ad))
        for ls in (KrylovJL_GMRES(), KLUFactorization(), nothing)
            alg = ls === nothing ? FBDF(nlsolve = nsa) : FBDF(linsolve = ls, nlsolve = nsa)
            sol = solve(prob, alg; reltol = 1.0e-5, abstol = 1.0e-8)
            @test SciMLBase.successful_retcode(sol)
            @test norm(sol.u[end] .- refsol.u[end]) / norm(refsol.u[end]) < 1.0e-5
        end
    end
end

@testset "Krylov: the ODE WOperator is reused matrix-free" begin
    integ = init(
        prob,
        FBDF(
            linsolve = KrylovJL_GMRES(),
            nlsolve = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff()))
        );
        reltol = 1.0e-5, abstol = 1.0e-8
    )
    c = integ.cache.nlsolver.cache
    @test c.W !== nothing
    # The inner NonlinearSolve holds the ODE's WOperator as its Jacobian (matrix-free),
    # not a residual-derived JacobianOperator.
    @test c.cache.jac_cache.J isa AbstractSciMLOperator
end
