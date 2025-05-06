using OrdinaryDiffEq, OrdinaryDiffEqCore, Test, LinearAlgebra
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear
using DiffEqDevTools, ADTypes

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
v = @inferred OrdinaryDiffEqCore.ode_interpolant(
    1.0, integrator1, integrator1.opts.save_idxs,
    Val{0})
@inferred OrdinaryDiffEqCore.ode_interpolant!(
    v, 1.0, integrator1, integrator1.opts.save_idxs,
    Val{0})
v = @inferred OrdinaryDiffEqCore.ode_extrapolant(
    1.0, integrator1, integrator1.opts.save_idxs,
    Val{0})
@inferred OrdinaryDiffEqCore.ode_extrapolant!(
    v, 1.0, integrator1, integrator1.opts.save_idxs,
    Val{0})

@testset "Mixed adaptivity" begin
    reverse_choice(integrator) = (Int(integrator.t > 0.5) + 1)
    alg_mixed = CompositeAlgorithm((Tsit5(), ABM54()), choice_function)
    alg_mixed_r = CompositeAlgorithm((ABM54(), Tsit5()), reverse_choice)
    alg_mixed2 = CompositeAlgorithm((Tsit5(), ABM54()), reverse_choice)

    @test_throws ArgumentError solve(prob_ode_linear, alg_mixed)
    sol2 = solve(prob_ode_linear, Tsit5())
    sol3 = solve(prob_ode_linear, alg_mixed; dt = 0.05, adaptive = false)
    sol4 = solve(prob_ode_linear, alg_mixed_r; dt = 0.05, adaptive = false)
    sol5 = solve(prob_ode_linear, alg_mixed2; dt = 0.05, adaptive = false)
    @test sol3.t == sol4.t && sol3.u == sol4.u
    @test sol3(0.8)≈sol2(0.8) atol=1e-4
    @test sol5(0.8)≈sol2(0.8) atol=1e-4
end

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

sol = solve(prob,
    alg = OrdinaryDiffEqCore.AutoAlgSwitch(ExplicitRK(constructVerner7()), Rodas5()))
@test sol.t[end] == 1000.0

prob = remake(prob_ode_2Dlinear, u0 = rand(ComplexF64, 2, 2))
sol = solve(prob, AutoTsit5(Rosenbrock23(autodiff = AutoFiniteDiff()))) # Complex and AD don't mix
@test sol.retcode == ReturnCode.Success

# https://github.com/SciML/ModelingToolkit.jl/issues/3043
function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    nothing
end
M = [1.0 0 0
     0 1.0 0
     0 0 0]
f = ODEFunction(rober, mass_matrix = M)
prob_mm = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1e5), (0.04, 3e7, 1e4))
cb = DiscreteCallback(
    (u, t, integrator) -> true, (integrator) -> u_modified!(integrator, true))
sol = solve(prob_mm, DefaultODEAlgorithm(), callback = cb)
