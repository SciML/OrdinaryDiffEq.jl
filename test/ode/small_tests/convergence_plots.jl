using Revise, OrdinaryDiffEq, DiffEqDevTools, Test, Random, ODEInterfaceDiffEq
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear
probArr = Vector{ODEProblem}(undef, 2)
plotArr = Vector{Plots.Plot{Plots.GRBackend}}(undef, 2)
probArr[1] = prob_ode_linear
probArr[2] = prob_ode_2Dlinear
titleArr = ["ExtrapolationMidpointDeuflhard: Linear ODE", "ExtrapolationMidpointDeuflhard: 2D Linear ODE"]
Random.seed!(100)
## Convergence Testing
 dts = 1 .//2 .^(12:-1:1) |>sort
 # dts = (100:10:1000)//1000 |>collect
testTol = 0.2
for i = 1:2
  global dts
  prob = probArr[i]
  plt = plot(legend = :topleft, title = titleArr[i], xlabel = "Stepsize", ylabel = "Final Error")
  # Order Convergence test
  alg2 = DP5()
  sim2 = test_convergence(dts,prob,alg2)
  plot!(plt,dts,markershape = :xcross, sim2.errors[:final], axis = :log10, label = "DP5")
  plotArr[i] = plt
  for j = 1:6, symbol in [:romberg]
      alg = ExtrapolationMidpointDeuflhard(min_extrapolation_order=j,
      init_extrapolation_order=j, max_extrapolation_order=j, sequence_symbol = symbol)
      sim = test_convergence(dts,prob,alg, dense = false)
      plot!(plt,dts, markershape = :circle,sim.errors[:final], axis = :log10, label = "Extr.Ord = $j")
  end
display(plt)
end
