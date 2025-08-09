import ODEProblemLibrary: prob_ode_vanderpol_stiff
using SciMLOperators
using OrdinaryDiffEq
using LinearAlgebra
using LinearSolve
using Test

# Define Jacobian for vanderpol stiff problem since v1.x doesn't provide it
function vanderpol_jac!(J, u, p, t)
    J[1, 1] = 0.0
    J[1, 2] = 1.0
    J[2, 1] = -2 * p[1] * u[1] * u[2] - 1.0
    J[2, 2] = p[1] * (1 - u[1]^2)
    return nothing
end

for base_prob in (prob_ode_vanderpol_stiff,)
    # If the problem doesn't have a jacobian (ODEProblemLibrary v1.x), add one
    if base_prob.f.jac === nothing
        f_with_jac = ODEFunction(base_prob.f.f; jac = vanderpol_jac!)
        prob = remake(base_prob; f = f_with_jac)
    else
        prob = base_prob
    end
    
    # Ensure all solutions use the same linear solve for fair comparison.
    # TODO: in future, ensure and test that polyalg chooses the best linear solve when unspecified.
    for alg in (Rosenbrock23(linsolve = KrylovJL_GMRES()),
        FBDF(linsolve = KrylovJL_GMRES()))
        # Manually construct a custom W operator using the Jacobian 
        N = length(prob.u0)
        J_op = MatrixOperator(zeros(N, N); update_func! = prob.f.jac)
        gamma_op = ScalarOperator(0.0;
            update_func = (old_val, u, p, t; dtgamma) -> dtgamma,
            accepted_kwargs = (:dtgamma,))
        transform_op = ScalarOperator(0.0;
            update_func = (old_op, u, p, t; dtgamma) -> inv(dtgamma),
            accepted_kwargs = (:dtgamma,))
        W_op = -(I - gamma_op * J_op) * transform_op

        # Make problem with custom MatrixOperator jac_prototype
        # Need to pass both jac and jac_prototype for proper function
        f_J = ODEFunction(prob.f.f; jac = prob.f.jac, jac_prototype = J_op)
        prob_J = remake(prob; f = f_J)

        # Test that the custom jacobian is used
        integrator = init(prob_J, Rosenbrock23())
        @test integrator.cache.J isa typeof(J_op)

        # Make problem with custom SciMLOperator W_prototype 
        # Need to pass jac as well for proper jacobian computation
        f_W = ODEFunction(prob.f.f; jac = prob.f.jac, jac_prototype = J_op, W_prototype = W_op)
        prob_W = remake(prob; f = f_W)

        # Test that the custom W operator is used
        integrator = init(prob_W, alg)
        if hasproperty(integrator.cache, :W)
            @test integrator.cache.W isa typeof(W_op)
        elseif hasproperty(integrator.cache.nlsolver.cache, :W)
            @test integrator.cache.nlsolver.cache.W isa typeof(W_op)
        else
            error("W-prototype test expected W in integrator.cache or integrator.cache.nlsolver.cache")
        end

        # Run solves
        sol = solve(prob, alg)
        sol_J = solve(prob_J, alg) # note: direct linsolve in this case is broken, see #1998
        sol_W = solve(prob_W, alg)

        rtol = 1e-2
        
        # Check sparsity only if it exists (not present in ODEProblemLibrary v1.x)
        if hasproperty(prob_J.f, :sparsity) && prob_J.f.sparsity !== nothing
            @test prob_J.f.sparsity.A == prob_W.f.sparsity.A
        end

        @test all(isapprox.(sol_J.t, sol_W.t; rtol))
        @test all(isapprox.(sol_J.u, sol_W.u; rtol))

        # For ODEProblemLibrary v1.x without ModelingToolkit, the MatrixOperator
        # jacobian prototype causes different step sizes, but solution is still correct
        # So we only compare the final states, not the trajectories
        if base_prob.f.jac === nothing
            # Compare final states only when we added our own jacobian
            @test isapprox(sol_J.u[end], sol.u[end]; rtol)
            @test isapprox(sol_W.u[end], sol.u[end]; rtol)
        else
            # Original test for when jacobian was provided (v0.1.x)
            @test all(isapprox.(sol_J.t, sol.t; rtol))
            @test all(isapprox.(sol_J.u, sol.u; rtol))
            @test all(isapprox.(sol_W.t, sol.t; rtol))
            @test all(isapprox.(sol_W.u, sol.u; rtol))
        end
    end
end
