using StochasticDiffEq, DiffEqDevTools, Test, Random
using SDEProblemLibrary: prob_sde_wave, prob_sde_cubic, prob_sde_additive
Random.seed!(100)
dts = (1 / 2) .^ (10:-1:2) #14->7 good plot

print("prob_sde_wave")
prob = prob_sde_wave
solve(prob, ImplicitEM(), dt = 0.1)
sim = test_convergence(dts, prob, ImplicitEM(), trajectories = Int(1.0e1))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2
sim = test_convergence(
    dts, prob, ImplicitEM(nlsolve = StochasticDiffEq.NLFunctional()), trajectories = Int(1.0e1)
)
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2
sim = test_convergence(dts, prob, ImplicitRKMil(), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 1) < 0.2
sim = test_convergence(dts, prob, EM(), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2
sim = test_convergence(dts, prob, ISSEM(), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2
sim = test_convergence(dts, prob, ISSEM(nlsolve = StochasticDiffEq.NLFunctional()), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2
sim = test_convergence(dts, prob, LambaEM(), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2
sim2 = test_convergence(dts, prob, RKMil(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, RKMilCommute(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, RKMilGeneral(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
print(".")
sim2 = test_convergence(dts, prob, WangLi3SMil_A(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_B(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_C(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_D(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_E(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_F(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
print(".")
eigen_est = (integrator) -> integrator.eigen_est = 10.0
sim2 = test_convergence(dts, prob, SROCK1(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SROCK1(eigen_est = eigen_est), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SROCK2(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SROCK2(eigen_est = eigen_est), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SROCKEM(strong_order_1 = false), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 0.5) < 0.2
sim2 = test_convergence(dts, prob, SROCKEM(eigen_est = eigen_est), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SROCKEM(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SKSROCK(eigen_est = eigen_est), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 0.5) < 0.2
sim2 = test_convergence(dts, prob, SKSROCK(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 0.5) < 0.2
sim2 = test_convergence(dts, prob, SROCKC2(eigen_est = eigen_est), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SROCKC2(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2

#omitting tests for incomplete methods
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2

sim3 = test_convergence(dts, prob, SRI(), trajectories = Int(1.0e1))
@test abs(sim3.𝒪est[:final] - 1.5) < 0.3
sim4 = test_convergence(dts, prob, SRIW1(), trajectories = Int(1.0e1))
@test abs(sim4.𝒪est[:final] - 1.5) < 0.3
sim5 = test_convergence(dts, prob, SRIW2(), trajectories = Int(1.0e1))
@test abs(sim5.𝒪est[:final] - 1.5) < 0.3
sim6 = test_convergence(dts, prob, SOSRI(), trajectories = Int(1.0e1))
@test abs(sim6.𝒪est[:final] - 1.5) < 0.3
sim7 = test_convergence(dts, prob, SOSRI2(), trajectories = Int(1.0e1))
@test abs(sim7.𝒪est[:final] - 1.5) < 0.3
println()

print("prob_sde_cubic")
prob = prob_sde_cubic
sim = test_convergence(dts, prob, EM(), trajectories = Int(1.0e1))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2
sim = test_convergence(dts, prob, LambaEM(), trajectories = Int(1.0e1))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2
sim = test_convergence(dts, prob, ISSEM(), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2
sim = test_convergence(dts, prob, ImplicitEM(), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2
sim = test_convergence(dts, prob, ImplicitRKMil(), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 1) < 0.2
sim2 = test_convergence(dts, prob, RKMil(), trajectories = Int(1.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.22
sim2 = test_convergence(dts, prob, RKMilCommute(), trajectories = Int(1.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.22
sim2 = test_convergence(dts, prob, RKMilGeneral(), trajectories = Int(1.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.22
print(".")
sim2 = test_convergence(dts, prob, WangLi3SMil_A(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_B(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_C(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_D(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_E(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_F(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
print(".")
sim2 = test_convergence(dts, prob, SROCK1(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SROCK2(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SROCKEM(strong_order_1 = false), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 0.5) < 0.2
sim2 = test_convergence(dts, prob, SROCKEM(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SKSROCK(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 0.5) < 0.2
sim2 = test_convergence(dts, prob, SROCKC2(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2

#omitting tests for incomplete methods
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2

sim3 = test_convergence(dts, prob, SRI(), trajectories = Int(1.0e1))
@test abs(sim3.𝒪est[:final] - 1.5) < 0.3
sim4 = test_convergence(dts, prob, SRIW1(), trajectories = Int(1.0e1))
@test abs(sim4.𝒪est[:final] - 1.5) < 0.3
sim5 = test_convergence(dts, prob, SRIW2(), trajectories = Int(1.0e1))
@test abs(sim5.𝒪est[:final] - 1.5) < 0.3
sim6 = test_convergence(dts, prob, SOSRI(), trajectories = Int(1.0e1))
@test abs(sim6.𝒪est[:final] - 1.5) < 0.3
sim7 = test_convergence(dts, prob, SOSRI2(), trajectories = Int(1.0e1))
@test abs(sim7.𝒪est[:final] - 1.5) < 0.3
println()

print("prob_sde_additive")
prob = prob_sde_additive
sim = test_convergence(dts, prob, EM(), trajectories = Int(1.0e1))
@test abs(sim.𝒪est[:l2] - 1) < 0.2
sim = test_convergence(dts, prob, LambaEM(), trajectories = Int(1.0e1))
@test abs(sim.𝒪est[:l2] - 1) < 0.2
dts = (1 / 2) .^ (10:-1:1)
sim = test_convergence(dts, prob, ISSEM(), trajectories = Int(1.0e3))
@test abs(sim.𝒪est[:l2] - 1) < 0.2
sim = test_convergence(dts, prob, ImplicitEM(), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 1) < 0.2
sim = test_convergence(dts, prob, ImplicitRKMil(), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 1) < 0.2
sim2 = test_convergence(dts, prob, RKMil(), trajectories = Int(1.0e1))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, RKMilCommute(), trajectories = Int(1.0e1))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, RKMilGeneral(), trajectories = Int(1.0e1))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
print(".")
sim2 = test_convergence(dts, prob, WangLi3SMil_A(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_B(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_C(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_D(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_E(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, WangLi3SMil_F(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
print(".")
sim2 = test_convergence(dts, prob, SROCK1(), trajectories = Int(1.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SROCK2(), trajectories = Int(1.0e2))
@test abs(sim2.𝒪est[:l∞] - 2) < 0.2
@time sim2 = test_convergence(dts, prob, SROCKEM(strong_order_1 = false), trajectories = Int(1.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
@time sim2 = test_convergence(dts, prob, SROCKEM(), trajectories = Int(1.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
sim2 = test_convergence(dts, prob, SKSROCK(), trajectories = Int(1.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2

dts = (1 / 2) .^ (10:-1:4)
sim2 = test_convergence(dts, prob, SROCKC2(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:l∞] - 1) < 0.2
dts = (1 / 2) .^ (10:-1:2)

#omitting tests for incomplete methods
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),trajectories=Int(2e2))
# @test abs(sim.𝒪est[:l∞]- 1) < 0.2

sim3 = test_convergence(dts, prob, SRI(), trajectories = Int(1.0e1))
@test abs(sim3.𝒪est[:final] - 2) < 0.3
sim3 = test_convergence(dts, prob, SRIW2(), trajectories = Int(1.0e1))
@test abs(sim3.𝒪est[:final] - 2) < 0.3

sim3 = test_convergence(dts, prob, SOSRI(), trajectories = Int(1.0e1))
@test abs(sim3.𝒪est[:final] - 2) < 0.3
sim3 = test_convergence(dts, prob, SOSRI2(), trajectories = Int(1.0e1))
@test abs(sim3.𝒪est[:final] - 2) < 0.3
print(".")
sim4 = test_convergence(dts, prob, SRA(), trajectories = Int(1.0e1))
@test abs(sim4.𝒪est[:final] - 2) < 0.3
sim5 = test_convergence(dts, prob, SRA1(), trajectories = Int(1.0e1))
@test abs(sim5.𝒪est[:final] - 2) < 0.3
sim6 = test_convergence(dts, prob, SRA2(), trajectories = Int(1.0e1))
@test abs(sim6.𝒪est[:final] - 2) < 0.3
sim7 = test_convergence(dts, prob, SRA3(), trajectories = Int(1.0e1))
@test abs(sim7.𝒪est[:final] - 2.5) < 0.3
sim8 = test_convergence(dts, prob, SOSRA(), trajectories = Int(1.0e1))
@test abs(sim8.𝒪est[:final] - 2) < 0.3
sim9 = test_convergence(dts, prob, SOSRA2(), trajectories = Int(1.0e1))
@test abs(sim9.𝒪est[:final] - 2) < 0.3
print(".")
dts = (1 / 2) .^ (10:-1:5) #14->7 good plot
sim10 = test_convergence(dts, prob, SKenCarp(), trajectories = Int(1.0e1))
@test abs(sim10.𝒪est[:final] - 2) < 0.3
sim10 = test_convergence(
    dts, prob, SKenCarp(nlsolve = StochasticDiffEq.NLFunctional()), trajectories = Int(1.0e1)
)
@test abs(sim10.𝒪est[:final] - 2) < 0.3

sim2 = test_convergence(dts, prob, SRA(tableau = StochasticDiffEq.constructSRA2()), trajectories = Int(1.0e1))
@test abs(sim2.𝒪est[:final] - 2) < 0.3
sim3 = test_convergence(dts, prob, SRA(tableau = StochasticDiffEq.constructSRA3()), trajectories = Int(1.0e2))
@test abs(sim3.𝒪est[:final] - 2.0) < 0.3
sim6 = test_convergence(dts, prob, SRIW1(), trajectories = Int(1.0e1))
@test abs(sim6.𝒪est[:final] - 2) < 0.3
sim2 = test_convergence(
    dts, prob, SRA(tableau = StochasticDiffEq.constructExplicitSKenCarp()),
    trajectories = Int(1.0e1)
)
@test abs(sim2.𝒪est[:final] - 2) < 0.3
println(".")
