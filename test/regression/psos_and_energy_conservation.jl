using OrdinaryDiffEq, ADTypes, Test, Random, LinearAlgebra, SparseArrays

# Parameters
Nc = 22
η = 1.0
κ = 1.0
T = (0.0, 100.0)

# Matrix definitions
A = sparse(diagm(1 => sqrt.(Complex[1:Nc;])))
H = η * (A + A') - 1.0im * κ * A' * A

u0 = zeros(ComplexF64, Nc + 1)
u0[1] = 1.0
function f_psos(du, u, t, p)
    du .= -1.0im * H * u
end

# Callback
rng = MersenneTwister(rand(UInt))
jumpnorm = Ref(rand(rng))
djumpnorm(x::Vector{ComplexF64}, t, integrator) = norm(x)^2 - (1 - jumpnorm[])
function dojump(integrator)
    x = integrator.u
    t = integrator.t

    x .= normalize(A * x)
    jumpnorm[] = rand(rng)
end

cb = ContinuousCallback(djumpnorm, dojump)

prob = ODEProblem{true}(f_psos, u0, T)

sol_tot = []
Ntraj = 100
for i in 1:Ntraj
    rng = MersenneTwister(rand(UInt))
    # Tweaking tolerances and dtmax also is not reliable
    sol = solve(prob, DP5(), save_everystep = true, callback = cb,
        abstol = 1e-8, reltol = 1e-6, dtmax = 10)
    push!(sol_tot, sol)
end

# This number has to be η^2/κ^2 in steady-state; all trajectories should converge there
n = [[norm(A * normalize(s.u[j]))^2 for j in 1:length(s.t)] for s in sol_tot]

@test all(η^2 / κ^2 .≈ [k[end] for k in n])

#=
using Plots
gr()

p1 = plot(sol_tot[1].t, n[1], lw = 2)
for i=2:Ntraj
    plot!(p1,sol_tot[i].t, n[i])
end

p2 = plot(sol_tot[1].t, norm.(sol_tot[1].u).^2)
for i=2:Ntraj
    plot!(p2,sol_tot[i].t, norm.(sol_tot[i].u).^2)
end
=#

using OrdinaryDiffEq, DiffEqCallbacks, Test

# Initial state
u0 = [0, -0.25, 0.42081, 0]

function hheom!(du, u, p, t)
    du[1] = u[3]
    du[2] = u[4]
    du[3] = -u[1] - 2u[1] * u[2]
    du[4] = -u[2] - (u[1]^2 - u[2]^2)
    return nothing
end

@inline Vhh(q1, q2) = 1 // 2 * (q1^2 + q2^2 + 2q1^2 * q2 - 2 // 3 * q2^3)
@inline Thh(p1, p2) = 1 // 2 * (p1^2 + p2^2)
@inline Hhh(q1, q2, p1, p2) = Thh(p1, p2) + Vhh(q1, q2)
@inline Hhh(u::AbstractVector) = Hhh(u...)

# Energy
const E = Hhh(u0)

function ghh(resid, u, p)
    resid[1] = Hhh(u[1], u[2], u[3], u[4]) - E
    resid[2:4] .= 0
end

# energy conserving callback:
# important to use save = false, I don't want rescaling points
cb = ManifoldProjection(ghh, abstol = 1e-13, save = false, autodiff = AutoForwardDiff())

# Callback for Poincare surface of section
function psos_callback(j, direction = +1, offset::Real = 0,
        callback_kwargs = Dict{Symbol, Any}(:abstol => 1e-9))

    # Prepare callback:
    s = sign(direction)
    cond = (u, t, integrator) -> s * (u - offset)
    affect! = (integrator) -> nothing

    cb = DiffEqBase.ContinuousCallback(cond, nothing, affect!; callback_kwargs...,
        save_positions = (true, false), idxs = j)
end

# with this callback, the saved values of variable 1 should be zero
poincarecb = psos_callback(1)

totalcb = CallbackSet(poincarecb, cb)

prob = ODEProblem(hheom!, u0, (0.0, 100.0), callback = totalcb)

extra_kw = Dict(:save_start => false, :save_end => false)
DEFAULT_DIFFEQ_KWARGS = Dict{Symbol, Any}(:abstol => 1e-9, :reltol => 1e-9)

sol = solve(prob, Vern9(); extra_kw..., DEFAULT_DIFFEQ_KWARGS..., save_everystep = false)

Es = [Hhh(sol[:, i]) for i in 1:length(sol)]
Eerror = maximum(@. abs(E - Es))

a = sol[1, :]

@test Eerror < 1e-10
for el in a
    @test abs(el) < 1e-10
end
