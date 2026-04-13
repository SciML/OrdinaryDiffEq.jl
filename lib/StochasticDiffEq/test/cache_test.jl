using StochasticDiffEq, Random

function f(du, u, p, t)
    for i in 1:length(u)
        du[i] = (0.2 / length(u)) * u[i]
    end
    return
end

function g(du, u, p, t)
    for i in 1:length(u)
        du[i] = (0.2 / length(u)) * u[i]
    end
    return
end

function condition(u, t, integrator)
    return 1 - maximum(u)
end

function callback_affect!(integrator)
    u = integrator.u
    resize!(integrator, length(u) + 1)
    maxidx = findmax(u)[2]
    Θ = rand()
    u[maxidx] = Θ
    u[end] = 1 - Θ
    return nothing
end

callback = ContinuousCallback(condition, callback_affect!)

u0 = [0.2]
tspan = (0.0, 100.0)
prob = SDEProblem(f, g, u0, tspan)

Random.seed!(3)
sol = solve(prob, SRIW1(), callback = callback)

#=
using Plots; pyplot()
p1 = plot(sol,vars=(0,1),plotdensity=10000,title="Amount of X in Cell 1")
scatter!(sol,denseplot=false)
p2 = plot(sol.t,map((x)->length(x),sol.u),lw=3,
     ylabel="Number of Cells",xlabel="Time")
plot(p1,p2,layout=(2,1),size=(600,1000))
=#

Random.seed!(3)
sol = solve(prob, EM(), callback = callback, dt = 1 / 4)
sol = solve(prob, RKMil(), callback = callback, dt = 1 / 4)

function g(du, u, p, t)
    for i in 1:length(u)
        du[i] = (0.3 / length(u))
    end
    return
end
prob = SDEProblem(f, g, u0, tspan)

sol = solve(prob, SRA1(), callback = callback)
