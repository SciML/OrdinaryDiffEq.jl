using OrdinaryDiffEq, DiffEqCallbacks, LinearAlgebra

# https://github.com/SciML/DiffEqBase.jl/issues/564 : Fixed
gravity = 9.8
stiffness = 500
equilibrium_length = 1
T = 5.0

f(u, p, t) = begin
    x1, x2, dx1, dx2 = u
    length = abs(x2 - x1)
    spring_force = stiffness * (equilibrium_length - length)
    ddx1 = -gravity - spring_force
    ddx2 = -gravity + spring_force
    if x1 <= 0
        ddx1 = max(0, ddx1)
    end
    if x2 <= 0
        ddx2 = max(0, ddx2)
    end
    [dx1, dx2, ddx1, ddx2]
end

sol = solve(
    ODEProblem(f, [5.0, 6.0, 0.0, 0.0], (0.0, T)),
    Rodas5P(),
    callback = ContinuousCallback(
        (u, _, _) -> u[1],
        (integrator) -> (integrator.u[1] = 0; integrator.u[3] = 0)
    ),
    # callback = ContinuousCallback((u, _, _) -> u[1], (integrator) -> (integrator.u[3] = 0)),
    reltol = 1.0e-5,
    abstol = 1.0e-5
)

@show sol.stats

# https://github.com/SciML/DiffEqBase.jl/issues/553 : Floating point issue is resolved but some other error occurs
function model(du, u, p, t)
    du[1] = 0.0
    for i in 2:(length(du) - 1)
        du[i] = p[i] * (u[i - 1] - u[i])
    end
    du[end] = p[end] * (p[1] * u[end - 1] - u[end])
    return nothing
end

perror = [
    1.0,
    0.02222434508140991,
    0.017030281542289794,
    0.015917011145559996,
    0.1608874463597176,
    0.13128016561792297,
    0.11056834258380167,
    0.5222141958458832,
    1.0711942201995688,
    0.2672878398678257,
    8.900058706990183,
    0.010760065201065117,
    0.016319181296867765,
    2.2693845639611925,
    0.2152216345154439,
    0.029186712540925457,
    0.21419429135100806,
    0.029177617589788596,
    0.03064986043089549,
    0.023280222517122397,
    6.931251277770224,
]
y_max = 0.002604806609572015
u0 = [1, zeros(length(perror) - 1)...]
tspan = (0.0, 5000.0)

condition(u, t, i) = (t == 1.0)
affect!(i) = (i.u[1] = 0.0)

condition2(u, t, i) = u[end] - y_max / 2.0
t_half_1 = 0.0
affect2!(i) = (t_half_1 = i.t)

prob = ODEProblem(model, u0, tspan, perror)
sol = solve(
    prob,
    Rosenbrock23();
    callback = CallbackSet(
        PositiveDomain(),
        DiscreteCallback(condition, affect!),
        ContinuousCallback(condition2, affect2!, terminate!)
    ),
    tstops = [1.0],
    force_dtmin = true
)

# https://github.com/SciML/DiffEqBase.jl/issues/515 : Fixed

using StaticArrays
using MultiScaleArrays

t_last = 0.0
function attactor(du, u, p, t)
    α, β = p
    n = length(u.nodes)
    return for k in 1:n
        du.nodes[k] = zero(du.nodes[k])
        for j in 1:n
            if (k == j)
                du.nodes[k] .+= [
                    u.nodes[k][3],
                    u.nodes[k][4],
                    -β * u.nodes[k][3],
                    -β * u.nodes[k][4],
                ]
            else
                du.nodes[k][3:4] .+= α * (u.nodes[j][1:2] - u.nodes[k][1:2])
            end
        end
    end
end

struct Thingy{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
end

struct PhysicsLaw{T <: AbstractMultiScaleArray, B <: Number} <:
    AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

Newton = construct(
    PhysicsLaw,
    [
        Thingy([-700.0, -350.0, 0.0, 0.0]),
        Thingy([-550.0, -150.0, 0.0, 0.0]),
        Thingy([-600.0, 15.0, 0.0, 10.0]),
        Thingy([200.0, -200.0, 5.0, -5.0]),
    ]
)

parameters = [1.0e-2, 0.06]

function condition(out, u, t, integrator)
    i = 0
    n = length(u.nodes)
    for k in 1:n
        for l in (k + 1):n
            i += 1
            out[i] = sum(abs2, u.nodes[k][1:2] .- u.nodes[l][1:2]) - 10000
        end
    end
    return
end

function collision_affect!(integrator, idx)
    i = 0
    u = integrator.u
    n = length(u.nodes)
    return for k in 1:n
        for l in (k + 1):n
            i += 1
            if idx == i
                x₁ = u.nodes[k][1:2]
                v₁ = u.nodes[k][3:4]
                x₂ = u.nodes[l][1:2]
                v₂ = u.nodes[l][3:4]
                # https://stackoverflow.com/a/35212639
                v₁ = (v₁ - 2 / (1 + 1) * (dot(v₁ - v₂, x₁ - x₂) / sum(abs2, x₁ - x₂) * (x₁ - x₂)))
                v₂ = -(v₂ - 2 / (1 + 1) * (dot(v₂ - v₁, x₂ - x₁) / sum(abs2, x₂ - x₁) * (x₂ - x₁)))

                println("Collision handled.")

                m = (x₁ + x₂) / 2

                u.nodes[k][3:4] .= v₁
                u.nodes[l][3:4] .= v₂

                set_u!(integrator, u)
                println(sqrt(sum(abs2, x₁ .- x₂)) - 100, ":", v₁ ./ v₂)
                println(
                    norm(v₁), ":", norm(v₂), ":", integrator.t, ":",
                    integrator.t - t_last
                )
                global t_last = integrator.t
                break
            end
        end
    end
end

cback = VectorContinuousCallback(
    condition,
    collision_affect!,
    (x -> Int(((x - 1) * x) / 2))(length(Newton.nodes))
)

problemp = ODEProblem(attactor, Newton, (0.0, Inf), parameters)

world = init(problemp, AutoTsit5(Rosenbrock23()); save_everystep = false, callback = cback)

dt = 0.2

for i in 1:1000
    step!(world, dt)
end

## https://github.com/SciML/OrdinaryDiffEq.jl/issues/1528

function f!(out, u, p, t)
    out[1] = 0
    out[2] = u[3]
    return out[3] = -1.0 * (u[2] - u[1])
end
u0 = [0, 0, 1.0]
function cond!(out, u, t, i)
    out[1] = u[3]
    return nothing
end
function terminate_affect!(int, idx)
    return terminate!(int)
end
cb = VectorContinuousCallback(cond!, terminate_affect!, nothing, 1)

u0 = [0.0, 0.0, 1.0]
prob = ODEProblem(f!, u0, (0.0, 10.0); callback = cb)
soln = solve(prob, Tsit5())
@test soln.t[end] ≈ 4.712347213360699 atol = 1e-4

odefun = ODEFunction((u, p, t) -> [u[2], u[2] - p]; mass_matrix = [1 0; 0 0])
callback = PresetTimeCallback(0.5, integ -> (integ.p = -integ.p))
prob = ODEProblem(odefun, [0.0, -1.0], (0.0, 1), 1; callback)
#test that reinit happens for both FSAL and non FSAL integrators
@testset "dae re-init" for alg in [FBDF(), Rodas5P()]
    sol = solve(prob, alg)
    # test that the callback flipping p caused u[2] to get flipped.
    first_t = findfirst(isequal(0.5), sol.t)
    @test sol.u[first_t][2] == -sol.u[first_t + 1][2]
end

daefun = DAEFunction((du, u, p, t) -> [du[1] - u[2], u[2] - p])
prob = DAEProblem(
    daefun, [0.0, 0.0], [0.0, -1.0], (0.0, 1), 1;
    differential_vars = [true, false], callback
)
sol = solve(prob, DFBDF())
# test that the callback flipping p caused u[2] to get flipped.
first_t = findfirst(isequal(0.5), sol.t)
@test sol.u[first_t][2] == -sol.u[first_t + 1][2]
