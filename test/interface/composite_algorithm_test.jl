using OrdinaryDiffEq, Test, LinearAlgebra
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear
prob = prob_ode_2Dlinear
choice_function(integrator) = (Int(integrator.t < 0.5) + 1)
alg_double = CompositeAlgorithm((Tsit5(), Tsit5()), choice_function)
alg_double2 = CompositeAlgorithm((Vern6(), Vern6()), choice_function)
alg_switch = CompositeAlgorithm((Tsit5(), Vern7()), choice_function)

@time sol1 = solve(prob_ode_linear, alg_double)
@time sol2 = solve(prob_ode_linear, Tsit5())
@test sol1.t == sol2.t
@test sol1(0.8) == sol2(0.8)

alg_double_erk = CompositeAlgorithm((ExplicitRK(), ExplicitRK()), choice_function)
@time sol1 = solve(prob_ode_linear, alg_double_erk)
@time sol2 = solve(prob_ode_linear, ExplicitRK())
@test sol1.t == sol2.t
@test sol1(0.8) == sol2(0.8)

integrator1 = init(prob, alg_double2)
integrator2 = init(prob, Vern6())
solve!(integrator1)
solve!(integrator2)

@test integrator1.sol.t == integrator2.sol.t

sol = solve(prob, alg_switch)
@inferred DiffEqBase.__init(prob, alg_switch)
v = @inferred OrdinaryDiffEq.ode_interpolant(1.0, integrator1, integrator1.opts.save_idxs,
                                             Val{0})
@inferred OrdinaryDiffEq.ode_interpolant!(v, 1.0, integrator1, integrator1.opts.save_idxs,
                                          Val{0})
v = @inferred OrdinaryDiffEq.ode_extrapolant(1.0, integrator1, integrator1.opts.save_idxs,
                                             Val{0})
@inferred OrdinaryDiffEq.ode_extrapolant!(v, 1.0, integrator1, integrator1.opts.save_idxs,
                                          Val{0})

condition(u, t, integrator) = t == 192.0
function affect!(integrator)
    integrator.u[1] += 14000
    integrator.u[2] += 14000
end
A = [-0.027671669470584172 -0.0 -0.0 -0.0 -0.0 -0.0;
     -0.0 -0.05540281553537378 -0.0 -0.0 -0.0 -0.0;
     0.011534597161021629 0.011933539591245327 -0.24891886153387743 0.023054812171672122 0.0 0.0;
     0.0 0.0 0.17011732278405356 -0.023054812171672122 0.0 0.0;
     0.01613707230956254 0.04346927594412846 0.03148193084515083 0.0 -1.5621055510998967e9 7.293040577236404;
     0.0 0.0 0.0 0.0 1.559509231932001e9 -7.293040577236404]
prob = ODEProblem((du, u, p, t) -> mul!(du, A, u), zeros(6), (0.0, 1000), tstops = [192],
                  callback = DiscreteCallback(condition, affect!));
sol = solve(prob, alg = AutoVern7(Rodas5()))
@test sol.t[end] == 1000.0

sol = solve(prob, alg = OrdinaryDiffEq.AutoAlgSwitch(ExplicitRK(constructVerner7()), Rodas5()))
@test sol.t[end] == 1000.0