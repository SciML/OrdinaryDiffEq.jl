import ODEProblemLibrary: prob_ode_vanderpol_stiff
using SciMLOperators
using OrdinaryDiffEq
using LinearAlgebra
using Setfield # TODO: not a depenedency
using LinearSolve

prob = prob_ode_vanderpol_stiff

sol1 = solve(prob, Rosenbrock23())

# Manually construct a custom W operator using the Jacobian 
J_op = MatrixOperator(rand(2, 2); update_func=prob.f.jac)
gamma_op = ScalarOperator(0.0; update_func=(old_val, u, p, t; dtgamma) -> dtgamma, accepted_kwargs=(:dtgamma,))
# make transform op
transform_op = ScalarOperator(0.0; 
                              update_func = (old_op, u, p, t; dtgamma, transform) -> transform ? inv(dtgamma) : one(dtgamma),
                              accepted_kwargs=(:dtgamma, :transform))
W = -(I - gamma_op * J_op) * transform_op 
# TODO: need to concretize for fair comparison

begin
f2 = prob.f
f2 = @set! f2.W_prototype = W
f2 = @set! f2.jac_prototype = J_op # TODO: this is slow, even without setting W_prototype (maybe LinearSolve uses wrong alg?)
prob2 = remake(prob; f = f2)
end
sol2 = solve(prob2, Rosenbrock23())

@test length(sol1.u) == length(sol2.u)
@test sol1.u[end] ≈ sol2.u[end] 
