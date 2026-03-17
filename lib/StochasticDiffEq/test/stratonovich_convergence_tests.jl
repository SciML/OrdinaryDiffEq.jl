using StochasticDiffEq, Test, Random, DiffEqDevTools
using SDEProblemLibrary: prob_sde_linear_stratonovich, prob_sde_2Dlinear_stratonovich
import SciMLBase

Random.seed!(100)
dts = 1 ./ 2 .^ (10:-1:2) #14->7 good plot

prob = prob_sde_linear_stratonovich
sim = test_convergence(dts, prob, EulerHeun(), trajectories = Int(5.0e2))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(dts, prob, LambaEulerHeun(), trajectories = Int(5.0e2))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(dts, prob, ISSEulerHeun(), trajectories = Int(5.0e2))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(dts, prob, ImplicitEulerHeun(), trajectories = Int(5.0e2))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(dts, prob, ImplicitEulerHeun(theta = 1), trajectories = Int(5.0e2))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(dts, prob, ImplicitEulerHeun(symplectic = true), trajectories = Int(5.0e2))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(
    dts, prob, RKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    trajectories = Int(5.0e2)
)
@test abs(sim.𝒪est[:l2] - 1) < 0.2

sim = test_convergence(
    dts, prob,
    RKMilGeneral(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    trajectories = Int(5.0e2)
)
@test abs(sim.𝒪est[:l2] - 1) < 0.2

sim = test_convergence(
    dts, prob,
    ImplicitRKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    trajectories = Int(5.0e2)
)
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(
    dts, prob, SROCK1(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    trajectories = Int(2.0e2)
)
@test abs(sim.𝒪est[:l2] - 1) < 0.15

sim = test_convergence(dts, prob, KomBurSROCK2(), trajectories = Int(2.0e2))
@test abs(sim.𝒪est[:final] - 2) < 0.2

println("Now 2D")

prob = prob_sde_2Dlinear_stratonovich

sim = test_convergence(dts, prob, EulerHeun(), trajectories = Int(5.0e1))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(dts, prob, LambaEulerHeun(), trajectories = Int(5.0e1))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(dts, prob, ISSEulerHeun(), trajectories = Int(5.0e1))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(dts, prob, ImplicitEulerHeun(), trajectories = Int(5.0e1))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(dts, prob, ImplicitEulerHeun(theta = 1), trajectories = Int(5.0e1))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(dts, prob, ImplicitEulerHeun(symplectic = true), trajectories = Int(5.0e1))
@test abs(sim.𝒪est[:l2] - 1) < 0.1

println("RKMils")

sim = test_convergence(
    dts, prob, RKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    trajectories = Int(1.0e2)
)
@test abs(sim.𝒪est[:l2] - 1) < 0.2

sim = test_convergence(
    dts, prob,
    RKMilGeneral(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    trajectories = Int(1.0e2)
)
@test abs(sim.𝒪est[:l2] - 1) < 0.2

sim = test_convergence(
    dts, prob,
    ImplicitRKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    trajectories = Int(1.0e2)
)
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(
    dts,
    prob,
    ImplicitRKMil(theta = 1, interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    trajectories = Int(1.0e2)
)
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(
    dts,
    prob,
    ImplicitRKMil(symplectic = true, interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    trajectories = Int(1.0e2)
)
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(
    dts, prob, SROCK1(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
    trajectories = Int(1.0e2)
)
@test abs(sim.𝒪est[:l2] - 1) < 0.1

sim = test_convergence(dts, prob, KomBurSROCK2(), trajectories = Int(2.0e2))
@test abs(sim.𝒪est[:final] - 2) < 0.2

Random.seed!(200)
sol = solve(prob, EulerHeun(), dt = 1 / 4)
