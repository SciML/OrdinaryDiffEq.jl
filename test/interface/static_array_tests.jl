using StaticArrays, Test
using OrdinaryDiffEq, OrdinaryDiffEqCore, OrdinaryDiffEqNonlinearSolve
using RecursiveArrayTools, ADTypes

u0 = VectorOfArray([fill(2, MVector{2, Float64}), ones(MVector{2, Float64})])
g0(u, p, t) = SA[u[1] + u[2], u[1]]
f = (du, u, p, t) -> begin
    for i in 1:2
        du[:, i] = g0(u[:, i], p, t)
    end
end
ode = ODEProblem(f, u0, (0.0, 1.0))
sol = solve(ode, Euler(), dt = 1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, Tsit5())
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, Vern9())
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)

u0 = VectorOfArray([fill(2, SVector{2, Float64}), ones(SVector{2, Float64})])
ode = ODEProblem(f, u0, (0.0, 1.0))
sol = solve(ode, Euler(), dt = 1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, Tsit5())
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, SSPRK22(), dt = 1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, ROCK4())
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)

u0 = ones(MVector{2, Float64})
ode = ODEProblem(g0, u0, (0.0, 1.0))
sol = solve(ode, Euler(), dt = 1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, Tsit5(), dt = 1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)

u0 = ones(SVector{2, Float64})
f = (u, p, t) -> u
ode = ODEProblem(f, u0, (0.0, 1.0))
sol = solve(ode, Euler(), dt = 1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, ImplicitEuler())
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, ImplicitEuler(nlsolve = OrdinaryDiffEqNonlinearSolve.NLAnderson()))
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, Tsit5(), dt = 1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)

#https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/373
function lorenz_static(u, p, t)
    dx = 10.0 * (u[2] - u[1])
    dy = u[1] * (28.0 - u[3]) - u[2]
    dz = u[1] * u[2] - (8 / 3) * u[3]
    @SVector [dx, dy, dz]
end

u0 = @SVector [1.0, 0.0, 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz_static, u0, tspan)
solve(prob, dt = 0.1, Rosenbrock23(autodiff = AutoFiniteDiff()))

# Check that ArrayPartitions of static vectors work
#https://github.com/SciML/OrdinaryDiffEq.jl/issues/1308
function lorenz_static(u::ArrayPartition, p, t)
    dx = 10.0 * (u[2] - u[1])
    dy = u[1] * (28.0 - u[3]) - u[2]
    dz = u[1] * u[2] - (8 / 3) * u[3]
    du1 = @SVector [dx, dy]
    du2 = @SVector [dz]
    ArrayPartition(du1, du2)
end

u01 = @SVector [1.0, 0.0]
u02 = @SVector [0.0]
u0ap = ArrayPartition(u01, u02)
probap = ODEProblem(lorenz_static, u0ap, tspan)

sol = solve(prob, dt = 1e-2, Heun())
solap = solve(probap, dt = 1e-2, Heun())
@test sol(30)≈solap(30) atol=1e-12

sol = solve(prob, dt = 1e-2, Tsit5())
solap = solve(probap, dt = 1e-2, Tsit5())
@test sol(30)≈solap(30) atol=1e-5

function rober(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    dy₁ = -k₁ * y₁ + k₃ * y₂ * y₃
    dy₂ = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    dy₃ = k₂ * y₂^2
    SA[dy₁, dy₂, dy₃]
end
prob = ODEProblem{false}(rober, SA[1.0, 0.0, 0.0], (0.0, 1e5), SA[0.04, 3e7, 1e4])
# Defaults to reltol=1e-3, abstol=1e-6
@test_nowarn sol = solve(
    prob, Rosenbrock23(autodiff = AutoForwardDiff(chunksize = 3)), save_everystep = false)
@test_nowarn sol = solve(
    prob, Rodas4(autodiff = AutoForwardDiff(chunksize = 3)), save_everystep = false)

function hires_4(u, p, t)
    y1, y2, y3, y4 = u
    dy1 = -1.71 * y1 + 0.43 * y2 + 8.32 * y3 + 0.0007
    dy2 = 1.71 * y1 - 8.75 * y2
    dy3 = -10.03 * y3 + 0.43 * y4 + 0.035 * y2
    dy4 = 8.32 * y2 + 1.71 * y3 - 1.12 * y4
    SA[dy1, dy2, dy3, dy4]
end

u0 = SA[1, 0, 0, 0.0057]
prob = ODEProblem(hires_4, u0, (0.0, 321.8122))
# Defaults to reltol=1e-3, abstol=1e-6
@test_nowarn sol = solve(
    prob, Rosenbrock23(autodiff = AutoForwardDiff(chunksize = 4)), save_everystep = false)
@test_nowarn sol = solve(
    prob, Rodas5(autodiff = AutoForwardDiff(chunksize = 4)), save_everystep = false)

function hires_5(u, p, t)
    y1, y2, y3, y4, y5 = u
    dy1 = -1.71 * y1 + 0.43 * y2 + 8.32 * y3 + 0.0007
    dy2 = 1.71 * y1 - 8.75 * y2
    dy3 = -10.03 * y3 + 0.43 * y4 + 0.035 * y5
    dy4 = 8.32 * y2 + 1.71 * y3 - 1.12 * y4
    dy5 = -1.745 * y5 + 0.43 * y2 + 0.43 * y4
    SA[dy1, dy2, dy3, dy4, dy5]
end

u0 = SA[1, 0, 0, 0, 0.0057]
prob = ODEProblem(hires_5, u0, (0.0, 321.8122))
# Defaults to reltol=1e-3, abstol=1e-6
@test_nowarn sol = solve(
    prob, Rosenbrock23(autodiff = AutoForwardDiff(chunksize = 5)), save_everystep = false)
@test_nowarn sol = solve(
    prob, Rodas4(autodiff = AutoForwardDiff(chunksize = 5)), save_everystep = false)

function hires(u, p, t)
    y1, y2, y3, y4, y5, y6, y7, y8 = u
    dy1 = -1.71 * y1 + 0.43 * y2 + 8.32 * y3 + 0.0007
    dy2 = 1.71 * y1 - 8.75 * y2
    dy3 = -10.03 * y3 + 0.43 * y4 + 0.035 * y5
    dy4 = 8.32 * y2 + 1.71 * y3 - 1.12 * y4
    dy5 = -1.745 * y5 + 0.43 * y6 + 0.43 * y7
    dy6 = -280.0 * y6 * y8 + 0.69 * y4 + 1.71 * y5 -
          0.43 * y6 + 0.69 * y7
    dy7 = 280.0 * y6 * y8 - 1.81 * y7
    dy8 = -280.0 * y6 * y8 + 1.81 * y7
    SA[dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8]
end

u0 = SA[1, 0, 0, 0, 0, 0, 0, 0.0057]
prob = ODEProblem(hires, u0, (0.0, 321.8122))
# Defaults to reltol=1e-3, abstol=1e-6
@test_nowarn sol = solve(
    prob, Rosenbrock23(autodiff = AutoForwardDiff(chunksize = 8)), save_everystep = false)
@test_nowarn sol = solve(
    prob, Rodas5(autodiff = AutoForwardDiff(chunksize = 8)), save_everystep = false)

const k1 = 0.35e0
const k2 = 0.266e2
const k3 = 0.123e5
const k4 = 0.86e-3
const k5 = 0.82e-3
const k6 = 0.15e5
const k7 = 0.13e-3
const k8 = 0.24e5
const k9 = 0.165e5
const k10 = 0.9e4
const k11 = 0.22e-1
const k12 = 0.12e5
const k13 = 0.188e1
const k14 = 0.163e5
const k15 = 0.48e7
const k16 = 0.35e-3
const k17 = 0.175e-1
const k18 = 0.1e9
const k19 = 0.444e12
const k20 = 0.124e4
const k21 = 0.21e1
const k22 = 0.578e1
const k23 = 0.474e-1
const k24 = 0.178e4
const k25 = 0.312e1

function pollu(y, p, t)
    r1 = k1 * y[1]
    r2 = k2 * y[2] * y[4]
    r3 = k3 * y[5] * y[2]
    r4 = k4 * y[7]
    r5 = k5 * y[7]
    r6 = k6 * y[7] * y[6]
    r7 = k7 * y[9]
    r8 = k8 * y[9] * y[6]
    r9 = k9 * y[11] * y[2]
    r10 = k10 * y[11] * y[1]
    r11 = k11 * y[13]
    r12 = k12 * y[10] * y[2]
    r13 = k13 * y[14]
    r14 = k14 * y[1] * y[6]
    r15 = k15 * y[3]
    r16 = k16 * y[4]
    r17 = k17 * y[4]
    r18 = k18 * y[16]
    r19 = k19 * y[16]
    r20 = k20 * y[17] * y[6]
    r21 = k21 * y[19]
    r22 = k22 * y[19]
    r23 = k23 * y[1] * y[4]
    r24 = k24 * y[19] * y[1]
    r25 = k25 * y[20]

    dy1 = -r1 - r10 - r14 - r23 - r24 +
          r2 + r3 + r9 + r11 + r12 + r22 + r25
    dy2 = -r2 - r3 - r9 - r12 + r1 + r21
    dy3 = -r15 + r1 + r17 + r19 + r22
    dy4 = -r2 - r16 - r17 - r23 + r15
    dy5 = -r3 + r4 + r4 + r6 + r7 + r13 + r20
    dy6 = -r6 - r8 - r14 - r20 + r3 + r18 + r18
    dy7 = -r4 - r5 - r6 + r13
    dy8 = r4 + r5 + r6 + r7
    dy9 = -r7 - r8
    dy10 = -r12 + r7 + r9
    dy11 = -r9 - r10 + r8 + r11
    dy12 = r9
    dy13 = -r11 + r10
    dy14 = -r13 + r12
    dy15 = r14
    dy16 = -r18 - r19 + r16
    dy17 = -r20
    dy18 = r20
    dy19 = -r21 - r22 - r24 + r23 + r25
    dy20 = -r25 + r24
    SA[dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8, dy9, dy10, dy11, dy12, dy13, dy14, dy15,
        dy16, dy17, dy18, dy19, dy20]
end

u0 = zeros(20)
u0[2] = 0.2
u0[4] = 0.04
u0[7] = 0.1
u0[8] = 0.3
u0[9] = 0.01
u0[17] = 0.007
u0 = SA[u0...]
prob = ODEProblem(pollu, u0, (0.0, 60.0))
@test_nowarn sol = solve(
    prob, Rosenbrock23(autodiff = AutoForwardDiff(chunksize = 8)), save_everystep = false)
@test_nowarn sol = solve(
    prob, Rodas5(autodiff = AutoForwardDiff(chunksize = 8)), save_everystep = false)

# DFBDF
g1(du, u, p, t) = du .^ 2 - conj.(u)
u0 = SA[-0.5048235596641171 - 0.8807809019469485im,
    -0.5086319184891589 - 0.8791877406854778im,
    -0.5015635095728721 + 0.8770989403497113im,
    -0.5031603140023296 - 0.8808350427797037im,
    -0.5033025831496138 + 0.877646898204028im,
    -0.5051689310722398 - 0.8753325698472263im,
    -0.4934794086951181 - 0.8648817236591336im,
    -0.5076312332711304 + 0.8736927915035966im,
    -0.4972563077522698 + 0.867122302026744im,
    -0.507021559962068 + 0.8853510586397142im]
du0 = SA[-0.5051593302918506 - 0.87178524227302im,
    -0.5035292188986685 - 0.8730255395885403im,
    -0.5043891667056177 + 0.8694664691998534im,
    -0.505596728341443 - 0.871084591593664im,
    -0.5041911921742974 + 0.8703512747251844im,
    -0.5027304380940516 - 0.8705784424498901im,
    -0.5011400789526957 - 0.8629141251757512im,
    -0.5014121747348508 + 0.8712321173163112im,
    -0.5011616766671037 + 0.8651123244481334im,
    -0.5065728050401669 + 0.8738635859036186im]
prob = DAEProblem(g1, du0, u0, (0.0, 10.0))
sol1 = solve(prob, DFBDF(autodiff = AutoFiniteDiff()), reltol = 1e-8, abstol = 1e-8)

g2(resid, du, u, p, t) = resid .= du .^ 2 - conj.(u)
prob = DAEProblem(g2, Array(du0), Array(u0), (0.0, 10.0))
sol2 = solve(prob, DFBDF(autodiff = AutoFiniteDiff()), reltol = 1e-8, abstol = 1e-8)

@test all(iszero, sol1[:, 1] - sol2[:, 1])
@test all(abs.(sol1[:, end] .- sol2[:, end]) .< 1.5e-6)
