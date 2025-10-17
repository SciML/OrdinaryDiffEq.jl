using OrdinaryDiffEq
using OrdinaryDiffEqFIRK
using ADTypes
using SparseConnectivityTracer
using LinearSolve

function test_sparse(ode_solver)
    function f(du, u, p, t)
        du .= [u[1], u[2]]
    end

    u0 = [1.0, 2.0]
    p = ()
    du0 = similar(u0)
    jac_prototype = float.(ADTypes.jacobian_sparsity(
        (du, u) -> f(du, u, p, 0.0),
        du0,
        u0, TracerSparsityDetector()))

    ode_fun = ODEFunction(f, jac_prototype=jac_prototype)
    prob = ODEProblem(ode_fun, u0, (0, 10))
    sol = solve(prob, ode_solver)
    return sol
end

println("Testing RadauIIA5 with KLUFactorization...")
test_sparse(RadauIIA5(;linsolve=KLUFactorization()))
println("Success!")
