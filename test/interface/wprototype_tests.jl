import ODEProblemLibrary: prob_ode_vanderpol_stiff
using SciMLOperators
using OrdinaryDiffEq
using LinearAlgebra
using LinearSolve
using Test

for prob in (prob_ode_vanderpol_stiff,)
    # Ensure all solutions use the same linear solve for fair comparison.
    # TODO: in future, ensure and test that polyalg chooses the best linear solve when unspecified.
    for alg in (
            Rosenbrock23(linsolve = KrylovJL_GMRES()),
            FBDF(linsolve = KrylovJL_GMRES()),
        )
        # Manually construct a custom W operator using the Jacobian
        N = length(prob.u0)
        J_op = MatrixOperator(zeros(N, N); update_func! = prob.f.jac)
        gamma_op = ScalarOperator(
            0.0;
            update_func = (old_val, u, p, t; gamma) -> gamma,
            accepted_kwargs = Val((:gamma,))
        )
        transform_op = ScalarOperator(
            0.0;
            update_func = (old_op, u, p, t; gamma) -> inv(gamma),
            accepted_kwargs = Val((:gamma,))
        )
        W_op = -(I - gamma_op * J_op) * transform_op

        # Make problem with custom MatrixOperator jac_prototype
        f_J = ODEFunction(prob.f.f; jac_prototype = J_op)
        prob_J = remake(prob; f = f_J)

        # Test that the custom jacobian is used
        integrator = init(prob_J, Rosenbrock23())
        @test integrator.cache.J isa typeof(J_op)

        # Make problem with custom SciMLOperator W_prototype
        f_W = ODEFunction(prob.f.f; jac_prototype = J_op, W_prototype = W_op)
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

        rtol = 1.0e-2

        @test prob_J.f.sparsity.A == prob_W.f.sparsity.A

        @test all(isapprox.(sol_J.t, sol_W.t; rtol))
        @test all(isapprox.(sol_J.u, sol_W.u; rtol))

        # Compare endpoint values. Different Jacobian computation methods (analytical
        # vs auto-diff) may cause slightly different adaptive step counts, so we
        # compare solution values rather than requiring identical time arrays.
        @test isapprox(sol_J.u[end], sol.u[end]; rtol)
        @test isapprox(sol_W.u[end], sol.u[end]; rtol)
    end
end
