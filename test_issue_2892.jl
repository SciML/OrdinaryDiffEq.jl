using OrdinaryDiffEq
using OrdinaryDiffEqFIRK
using ADTypes
using SparseConnectivityTracer
using LinearSolve
using Sparspak

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

println("Testing AdaptiveRadau with LUFactorization...")
test_sparse(AdaptiveRadau(;linsolve=LUFactorization()))       # Success
println("Success!")

println("Testing AdaptiveRadau with SparspakFactorization...")
test_sparse(AdaptiveRadau(;linsolve=SparspakFactorization())) # Success
println("Success!")

println("Testing QNDF with KLUFactorization...")
test_sparse(QNDF(;linsolve=KLUFactorization()))               # Success
println("Success!")

println("Testing AdaptiveRadau with KLUFactorization...")
test_sparse(AdaptiveRadau(;linsolve=KLUFactorization()))      # Error
println("Success!")
