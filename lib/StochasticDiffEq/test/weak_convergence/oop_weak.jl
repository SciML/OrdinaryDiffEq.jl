using StochasticDiffEq, DiffEqDevTools, Test
using Random
using SDEProblemLibrary: prob_sde_linear
Random.seed!(100)
dts = 1 .// 2 .^ (7:-1:3) #14->7 good plot

prob = prob_sde_linear
println("EM")
@time sim = test_convergence(
    dts, prob, EM(), save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("SimplifiedEM")
@time sim = test_convergence(
    dts, prob, SimplifiedEM(), save_everystep = false, trajectories = Int(1.0e6)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.35
println("RKMil")
sim = test_convergence(dts, prob, RKMil(), save_everystep = false, trajectories = Int(1.0e4))
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("RKMilCommute")
sim = test_convergence(
    dts, prob, RKMilCommute(), save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("RKMilGeneral")
sim = test_convergence(
    dts, prob, RKMilGeneral(), save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("SROCK1")
sim = test_convergence(dts, prob, SROCK1(), save_everystep = false, trajectories = Int(1.0e4))
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("SROCK2")
dts = 1 .// 2 .^ (8:-1:1) #14->7 good plot
sim = test_convergence(dts, prob, SROCK2(), save_everystep = false, trajectories = Int(7.0e5))
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35
#@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3
dts = 1 .// 2 .^ (7:-1:2) #14->7 good plot
println("SROCKEM")
sim = test_convergence(
    dts, prob, SROCKEM(strong_order_1 = false),
    save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("SROCKEM")
sim = test_convergence(
    dts, prob, SROCKEM(), save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("SKSROCK")
sim = test_convergence(
    dts, prob, SKSROCK(), save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3

#omitting tests for incomplete methods
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=1),save_everystep=false,trajectories=Int(4e4),
#                           weak_timeseries_errors=true)
# @test abs(sim.𝒪est[:weak_final]-2) < 0.3
# #@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
# #@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=2),save_everystep=false,trajectories=Int(4e4),
#                           weak_timeseries_errors=true)
# @test abs(sim.𝒪est[:weak_final]-2) < 0.3
# #@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
# #@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=3),save_everystep=false,trajectories=Int(4e4),
#                           weak_timeseries_errors=true)
# @test abs(sim.𝒪est[:weak_final]-2) < 0.3
# #@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
# #@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=4),save_everystep=false,trajectories=Int(4e4),
#                           weak_timeseries_errors=true)
# @test abs(sim.𝒪est[:weak_final]-2) < 0.3
# #@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
# #@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3
# sim = test_convergence(dts,prob,TangXiaoSROCK2(version_num=5),save_everystep=false,trajectories=Int(4e4),
#                           weak_timeseries_errors=true)
# @test abs(sim.𝒪est[:weak_final]-2) < 0.3
# #@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
# #@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3

println("WangLi3SMil_A")
dts = 1 .// 2 .^ (7:-1:3)
sim = test_convergence(
    dts, prob, WangLi3SMil_A(), save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("WangLi3SMil_B")
sim = test_convergence(
    dts, prob, WangLi3SMil_B(), save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("WangLi3SMil_C")
sim = test_convergence(
    dts, prob, WangLi3SMil_C(), save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("WangLi3SMil_D")
sim = test_convergence(
    dts, prob, WangLi3SMil_D(), save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("WangLi3SMil_E")
sim = test_convergence(
    dts, prob, WangLi3SMil_E(), save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("WangLi3SMil_F")
sim = test_convergence(
    dts, prob, WangLi3SMil_F(), save_everystep = false, trajectories = Int(1.0e4)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-1) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-1) < 0.3
println("SRI")
sim = test_convergence(dts, prob, SRI(), save_everystep = false, trajectories = Int(2.0e5))
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3
println("SRIW1")
sim = test_convergence(dts, prob, SRIW1(), save_everystep = false, trajectories = Int(1.0e5))
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
#@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3

#sim = test_convergence(dts,prob,DRI1(),save_everystep=false,trajectories=Int(4e4),
#                        weak_timeseries_errors=true)
#@test abs(sim.𝒪est[:weak_final]-2) < 0.3
##@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
##@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3
#sim = test_convergence(dts,prob,RI1(),save_everystep=false,trajectories=Int(4e4),
#                        weak_timeseries_errors=true)
#@test abs(sim.𝒪est[:weak_final]-2) < 0.3
##@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
##@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3
