using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using NonlinearSolve, OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, DiffEqDevTools
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using Test

@parameters k₁ k₂ k₃
@variables y₁(t) y₂(t) y₃(t)

eqs = [D(y₁) ~ -k₁ * y₁ + k₃ * y₂ * y₃,
    D(y₂) ~ k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃,
    D(y₃) ~ k₂ * y₂^2]
@mtkcompile rober = ODESystem(eqs, t)
prob = ODEProblem(rober, [[y₁, y₂, y₃] .=> [1.0; 0.0; 0.0]; [k₁, k₂, k₃] .=> (0.04, 3e7, 1e4)], (0.0, 1e5), jac = true)
prob2 = ODEProblem(rober, [[y₁, y₂, y₃] .=> [1.0; 0.0; 0.0]; [k₁, k₂, k₃] .=> (0.04, 3e7, 1e4)], (0.0, 1e5), jac = true, nlstep = true)

@test prob2.f.nlstep_data !== nothing

nlalg = NonlinearSolveAlg(NewtonRaphson(autodiff = AutoFiniteDiff()));
nlalgrobust = NonlinearSolveAlg(TrustRegion(autodiff = AutoFiniteDiff()));
sol1 = solve(prob, FBDF(autodiff=AutoFiniteDiff(), nlsolve = nlalg));
sol2 = solve(prob2, FBDF(autodiff=AutoFiniteDiff(), nlsolve = nlalg));

@test sol1.t != sol2.t
@test sol1 != sol2
@test sol1(sol1.t) ≈ sol2(sol1.t) atol=1e-3

sol1 = solve(prob, TRBDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalg));
sol2 = solve(prob2, TRBDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalg));

@test sol1.t != sol2.t
@test sol1.u != sol2.u
@test sol1(sol1.t) ≈ sol2(sol1.t) atol=1e-3

testprob = ODEProblem(rober, [[y₁, y₂, y₃] .=> [1.0; 0.0; 0.0]; [k₁, k₂, k₃] .=> (0.04, 3e7, 1e4)], (0.0, 1.0), nlstep = true)
@test testprob.f.nlstep_data !== nothing
sol2 = solve(testprob, TRBDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalg), adaptive=false, dt= 2.0^-15);

test_setup = Dict(:alg => FBDF(), :reltol => 1e-14, :abstol => 1e-14)
dts = 2.0 .^ (-10:-1:-15)
sim = analyticless_test_convergence(dts, testprob, TRBDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalgrobust), test_setup);
@test abs(sim.𝒪est[:l∞] - 2) < 0.2

dts = 2.0 .^ (-10:-1:-12)
sim = analyticless_test_convergence(dts, testprob, KenCarp4(autodiff=AutoFiniteDiff(), nlsolve = nlalgrobust), test_setup);
@test abs(sim.𝒪est[:l∞] - 4) < 0.2

dts = 2.0 .^ (-12:-1:-15)
sim = analyticless_test_convergence(dts, testprob, ABDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalgrobust), test_setup);
@test abs(sim.𝒪est[:l∞] - 2) < 0.2

dts = 2.0 .^ (-13:-1:-16)
sim = analyticless_test_convergence(dts, testprob, QNDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalgrobust), test_setup);
@test abs(sim.𝒪est[:l∞] - 2.5) < 0.2 # Superconvergence

dts = 2.0 .^ (-15:-1:-18)
sim = analyticless_test_convergence(dts, testprob, FBDF(autodiff=AutoFiniteDiff(), nlsolve = nlalgrobust), test_setup);
@test abs(sim.𝒪est[:l∞] - 1) < 0.3 # Only first order because adaptive order starts with Euler!

eqs_nonaut = [D(y₁) ~ -k₁ * y₁ + (t+1) * k₃ * y₂ * y₃,
    D(y₂) ~ k₁ * y₁ - (t+1) * k₂ * y₂^2 - (t+1) * k₃ * y₂ * y₃,
    D(y₃) ~ (t+1) * k₂ * y₂^2]
@mtkcompile rober_nonaut = ODESystem(eqs_nonaut, t)
prob = ODEProblem(rober_nonaut, [[y₁, y₂, y₃] .=> [1.0; 0.0; 0.0]; [k₁, k₂, k₃] .=> (0.04, 3e7, 1e4)], (0.0, 1e5), jac = true)
prob2 = ODEProblem(rober_nonaut, [[y₁, y₂, y₃] .=> [1.0; 0.0; 0.0]; [k₁, k₂, k₃] .=> (0.04, 3e7, 1e4)], (0.0, 1e5), jac = true, nlstep = true)

sol1 = solve(prob, FBDF(autodiff=AutoFiniteDiff(), nlsolve = nlalg));
sol2 = solve(prob2, FBDF(autodiff=AutoFiniteDiff(), nlsolve = nlalg));

@test sol1.t != sol2.t
@test sol1.u != sol2.u
@test sol1(sol1.t) ≈ sol2(sol1.t) atol=1e-3

sol1 = solve(prob, TRBDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalg));
sol2 = solve(prob2, TRBDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalg));

@test sol1.t != sol2.t
@test sol1 != sol2
@test sol1(sol1.t) ≈ sol2(sol1.t) atol=1e-4

testprob = ODEProblem(rober_nonaut, [[y₁, y₂, y₃] .=> [1.0; 0.0; 0.0]; [k₁, k₂, k₃] .=> (0.04, 3e7, 1e4)], (0.0, 1.0), nlstep = true)
@test testprob.f.nlstep_data !== nothing
sol2 = solve(testprob, TRBDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalg), adaptive=false, dt= 2.0^-15)

test_setup = Dict(:alg => FBDF(), :reltol => 1e-14, :abstol => 1e-14)
dts = 2.0 .^ (-10:-1:-15)
sim = analyticless_test_convergence(dts, testprob, TRBDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalgrobust), test_setup);
@test abs(sim.𝒪est[:l∞] - 2) < 0.2

dts = 2.0 .^ (-10:-1:-12)
sim = analyticless_test_convergence(dts, testprob, KenCarp4(autodiff=AutoFiniteDiff(), nlsolve = nlalgrobust), test_setup);
@test abs(sim.𝒪est[:l∞] - 4) < 0.2

dts = 2.0 .^ (-12:-1:-15)
sim = analyticless_test_convergence(dts, testprob, ABDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalgrobust), test_setup);
@test abs(sim.𝒪est[:l∞] - 2) < 0.2

dts = 2.0 .^ (-13:-1:-16)
sim = analyticless_test_convergence(dts, testprob, QNDF2(autodiff=AutoFiniteDiff(), nlsolve = nlalgrobust), test_setup);
@test abs(sim.𝒪est[:l∞] - 2.5) < 0.2 # Superconvergence

dts = 2.0 .^ (-15:-1:-18)
sim = analyticless_test_convergence(dts, testprob, FBDF(autodiff=AutoFiniteDiff(), nlsolve = nlalgrobust), test_setup);
@test abs(sim.𝒪est[:l∞] - 1) < 0.35 # Only first order because adaptive order starts with Euler!
