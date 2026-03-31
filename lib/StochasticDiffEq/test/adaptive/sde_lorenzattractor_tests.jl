using StochasticDiffEq, Random
using SDEProblemLibrary: prob_sde_lorenz
Random.seed!(100)
prob = prob_sde_lorenz

## Solve and plot
println("Plot 1")
sol1 = solve(prob, SRI(); dt = 1 / 2^(4), abstol = 10, reltol = 0)
#=
p1 = plot(sol1[:,1],sol1[:,2],sol1[:,3],title="Absolute Tolerance = 10",leg=false,
          top_margin=50px,right_margin=50px,xguide="X",yguide="Y",zguide="Z",bottom_margin=50px,
          guidefont=font(16),titlefont=font(18),tickfont=font(16))
=#

println("Plot 2")
sol2 = solve(prob, SRI(); dt = 1 / 2^(4), abstol = 1, reltol = 0)
#=
p2 = plot(sol2[:,1],sol2[:,2],sol2[:,3],title="Absolute Tolerance = 10",leg=false,
          top_margin=50px,right_margin=50px,xguide="X",yguide="Y",zguide="Z",bottom_margin=50px,
          guidefont=font(16),titlefont=font(18),tickfont=font(16))
=#

println("Plot 3")
sol3 = solve(prob, SRI(); dt = 1 / 2^(4), abstol = 1 / 10, reltol = 0)
#=
p3 = plot(sol3[:,1],sol3[:,2],sol3[:,3],title="Absolute Tolerance = 10",leg=false,
          top_margin=50px,right_margin=50px,xguide="X",yguide="Y",zguide="Z",bottom_margin=50px,
          guidefont=font(16),titlefont=font(18),tickfont=font(16))
=#

println("Plot 4")
sol4 = solve(prob, SRI(); dt = 1 / 2^(4), abstol = 1 / 100, reltol = 0)
#=
p4 = plot(sol4[:,1],sol4[:,2],sol4[:,3],title="Absolute Tolerance = 10",leg=false,
          top_margin=50px,right_margin=50px,xguide="X",yguide="Y",zguide="Z",bottom_margin=50px,
          guidefont=font(16),titlefont=font(18),tickfont=font(16))
=#

#plot(p1,p2,p3,p4,size=(1200,800),plot_title="Lorenz Attractors")
#gui()
