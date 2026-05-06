using DiffEqDevTools, Test, LinearAlgebra, DelayDiffEq
using OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqFIRK, OrdinaryDiffEqBDF
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error)) 
using BenchmarkTools

#tests against Hairer's RADAR problems
h(p, t; idxs = nothing) = 0.5

# state-dependent delay: τ(t) = y(t)
function delay(u, p, t)
    return u[1]
end

function f(du, u, h, p, t)
    τ = delay(u, p, t)
    du[1] = h(p, t - τ; idxs = 1)
end

# initial condition at t = 0
u0 = [1.0]
tspan = (0.0, 10.0)
p = nothing
prob = DDEProblem(f, u0, h, tspan, p; dependent_lags = (delay,))
sol = solve(prob, MethodOfSteps(Tsit5()))

# https://dieci.math.gatech.edu/preps/DieciLopez-Fili4.pdf   
# vector fields
# BELOW IS WRONG NEED TO FIX
function f1!(du, u, p, t)
    x1, x2 = u
    du[1] = x2
    du[2] = -x1 + 1/(1.2 - x2)
end

function f2!(du, u, p, t)
    x1, x2 = u
    du[1] = x2
    du[2] = -x1 - 1/(0.8 + x2)
end

# switching surface Σ: x2 = 0.2
condition(u, p, t) = u[2] - 0.2

# mode indicator (which vector field is active)
mode = Ref(1)

function f!(du, u, p, t)
    if mode[] == 1
        f1!(du, u, p, t)
    else
        f2!(du, u, p, t)
    end
end

# switch dynamics when crossing Σ
function affect!(integrator)
    mode[] = 2  # toggle 1 ↔ 2
end

cb = ContinuousCallback(condition, affect!, is_discontinuity = true;)
cb2 = ContinuousCallback(condition, affect!, is_discontinuity = false;)

u0 = [-0.4, -0.5]
tspan = (0.0, 10.0)

prob = ODEProblem(f!, u0, tspan, p)

sol_disco_tsit5 = solve(prob, Tsit5(), callback=cb)
sol_no_disco_tsit5 = solve(prob, Tsit5(), callback=cb2)

sol_disco_radau = solve(prob, RadauIIA5(), callback=cb)
sol_no_disco_radau = solve(prob, RadauIIA5(), callback=cb2)

sol_disco_rosenbrock = solve(prob, Rodas5P(), callback=cb)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(), callback=cb2)

sol_disco_bs3 = solve(prob, BS3(), callback=cb)
sol_no_disco_bs3 = solve(prob, BS3(), callback=cb2)
