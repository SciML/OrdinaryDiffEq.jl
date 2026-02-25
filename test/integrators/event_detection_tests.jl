using StaticArrays
using OrdinaryDiffEq
using DiffEqDevTools
using Test

@inbounds @inline function ż(z, p, t)
    A, B, D = p
    p₀, p₂ = z[1:2]
    q₀, q₂ = z[3:4]

    return SVector{4}(
        -A * q₀ - 3 * B / √2 * (q₂^2 - q₀^2) - D * q₀ * (q₀^2 + q₂^2),
        -q₂ * (A + 3 * √2 * B * q₀ + D * (q₀^2 + q₂^2)),
        A * p₀,
        A * p₂
    )
end

condition(u, t, integrator) = u
affect!(integrator) = nothing
function cbf(idx)
    return ContinuousCallback(
        condition,
        affect!, nothing, save_positions = (false, true), idxs = idx
    )
end
z0 = SVector{4}(7.1989885061904335, -0.165912283356219, 0.0, -3.63534900748947)

tspan = (0.0, 300.0)
prob = ODEProblem(ż, z0, tspan, (A = 1, B = 0.55, D = 0.4), callback = cbf(3))
sol = solve(
    prob, Vern9(), abstol = 1.0e-14, reltol = 1.0e-14,
    save_everystep = false, save_start = false, save_end = false, maxiters = 1.0e6
)

@test length(sol) > 100
@test SciMLBase.successful_retcode(sol)

prob = ODEProblem(ż, z0, (0, 400.0), (A = 1, B = 0.55, D = 0.4), callback = cbf(3))
sol = solve(
    prob, Vern9(), abstol = 1.0e-14, reltol = 1.0e-14, save_everystep = false,
    save_start = false, save_end = false, maxiters = 2.0e4
)

@test length(sol) > 100
@test SciMLBase.successful_retcode(sol)

prob = ODEProblem(ż, z0, (0, 5000.0), (A = 1, B = 0.55, D = 0.4), callback = cbf(3))
sol = solve(
    prob, Vern9(), abstol = 1.0e-14, reltol = 1.0e-14, save_everystep = false,
    save_start = false, save_end = false, maxiters = 1.0e6
)

@test length(sol) > 1500
@test SciMLBase.successful_retcode(sol)

@info "Bouncing Ball"

f = function (du, u, p, t)
    du[1] = u[2]
    return du[2] = -p[1]
end
function condition(u, t, integrator) # Event when event_f(u,t) == 0
    return u[1]
end
function affect!(integrator)
    return integrator.u[2] = -integrator.u[2]
end
cb2 = ContinuousCallback(condition, affect!)
tspan = (0.0, 10000.0)
u0 = [50.0, 0.0]
p = 9.8
prob = ODEProblem(f, u0, tspan, p)
sol = solve(prob, Tsit5(), callback = cb2)
@test minimum(Array(sol)) > -40

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/2055
for alg in (Rodas4(), Rodas4P(), Rodas5(), Rodas5P())
    sol2 = solve(prob, alg; callback = cb2)
    sol3 = appxtrue(sol, sol2)
    @test sol3.errors[:L2] < 1.0e-5
    @test sol3.errors[:L∞] < 5.0e-5
    @test sol3.errors[:final] < 2.0e-5
end

function fball(du, u, p, t)
    du[1] = u[2]
    du[2] = -p
    du[3] = u[4]
    return du[4] = 0.0
end
u0 = [50.0, 0.0, 0.0, 2.0]
tspan = (0.0, 15.0)
p = 9.8
prob = ODEProblem(fball, u0, tspan, p)

x = Ref(0)
y = Ref(0)
z = Ref(0)

function condition(out, u, t, integrator)
    out[1] = u[1]
    return out[2] = (10.0 - u[3])u[3]
end

function affect!(integrator, idx)
    return if idx == 1
        x[] += 1
        integrator.u[2] = -0.9integrator.u[2]
    elseif idx == 2
        y[] += 1
        integrator.u[4] = -0.9integrator.u[4]
    end
end

function affect_neg!(integrator, idx)
    z[] += 1
    @show integrator.u[1]
    return @show integrator.u[3]
end

cb = VectorContinuousCallback(condition, affect!, 2)
cb2 = VectorContinuousCallback(condition, affect_neg!, affect!, 2)

sol = solve(prob, Tsit5(), callback = cb, dt = 1.0e-3, adaptive = false)
@test x[] == 3
@test y[] == 2

sol2 = solve(prob, Tsit5(), callback = cb2, dt = 1.0e-3, adaptive = false)
@test sol.u == sol2.u
@test z[] == 0

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/1273

function du!(du, u, p, t)
    return du[1] = 1
end

callback = ContinuousCallback(
    (u, t, integrator) -> 1.0,
    (integrator) -> nothing
)

prob = ODEProblem(du!, [0], (0.0, 1.0), callback = callback)

solve(prob, Tsit5())
solve(prob, RadauIIA3())
solve(prob, RadauIIA5())

du(u, p, t) = [1.0]

prob = ODEProblem(du, [0], (0.0, 1.0), callback = callback)

solve(prob, Tsit5())
solve(prob, RadauIIA3())
solve(prob, RadauIIA5())
