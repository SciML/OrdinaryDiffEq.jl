using OrdinaryDiffEqSymplecticRK, Test, RecursiveArrayTools, DiffEqDevTools, Statistics
using OrdinaryDiffEqTsit5

u0 = fill(0.0, 2)
v0 = ones(2)
function f1_harmonic(dv, v, u, p, t)
    return dv .= -u
end
function f2_harmonic(du, v, u, p, t)
    return du .= v
end
function harmonic_analytic(y0, p, x)
    v0, u0 = y0.x
    return ArrayPartition(-u0 * sin(x) + v0 * cos(x), u0 * cos(x) + v0 * sin(x))
end
ff_harmonic = DynamicalODEFunction(f1_harmonic, f2_harmonic; analytic = harmonic_analytic)
prob = DynamicalODEProblem(ff_harmonic, v0, u0, (0.0, 5.0))

sol = solve(prob, SymplecticEuler(), dt = 1 / 2)
sol_verlet = solve(prob, VelocityVerlet(), dt = 1 / 100)
sol_ruth3 = solve(prob, Ruth3(), dt = 1 / 100)

interp_time = 0:0.001:5
interp = sol(0.5)
interps = sol(interp_time)

sol_tsit5 = solve(prob, Tsit5())

prob = SecondOrderODEProblem(f1_harmonic, v0, u0, (0.0, 5.0))

sol2 = solve(prob, SymplecticEuler(), dt = 1 / 2)
sol2_verlet = solve(prob, VelocityVerlet(), dt = 1 / 100)
sol2_ruth3 = solve(prob, Ruth3(), dt = 1 / 100)

sol2_verlet(0.1)

@test sol.u[end][1] == sol2.u[end][1]
@test sol_verlet.u[end][1] == sol2_verlet.u[end][1]
@test sol_ruth3.u[end][1] == sol2_ruth3.u[end][1]
@test sol.u[end][3] == sol2.u[end][3]
@test sol_verlet.u[end][3] == sol2_verlet.u[end][3]
@test sol_ruth3.u[end][3] == sol2_ruth3.u[end][3]

prob = DynamicalODEProblem(ff_harmonic, v0, u0, (0.0, 5.0))
println("Convergence tests")

dts = 1 .// 2 .^ (6:-1:3)
# Symplectic Euler
sim = test_convergence(dts, prob, SymplecticEuler(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 1 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 1 rtol = 1.0e-1
# Verlet
sim = test_convergence(dts, prob, VelocityVerlet(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1.0e-1
# Test that position converges faster for Verlet
position_error = :final => [
    mean(sim[i].u[2].x[1] - sim[i].u_analytic[2].x[1])
        for i in 1:length(sim)
]
@test first(DiffEqDevTools.calc𝒪estimates(position_error).second) ≈ 4.0 rtol = 1.0e-1

# 2nd Order Tableaus
sim = test_convergence(dts, prob, VerletLeapfrog(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1.0e-1
sim = test_convergence(dts, prob, LeapfrogDriftKickDrift(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1.0e-1
sim = test_convergence(dts, prob, PseudoVerletLeapfrog(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1.0e-1
sim = test_convergence(dts, prob, McAte2(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1.0e-1

# Ruth
sim = test_convergence(dts, prob, Ruth3(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 3 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 3 rtol = 1.0e-1
sim = test_convergence(dts, prob, McAte3(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 3 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 3 rtol = 1.0e-1
sim = test_convergence(dts, prob, CandyRoz4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, McAte4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, McAte42(), dense_errors = true)
@test_broken sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test_broken sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, CalvoSanz4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

dts = 1 .// 2 .^ (4:-1:0)
sim = test_convergence(dts, prob, McAte5(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

sim = test_convergence(dts, prob, Yoshida6(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4.69 rtol = 1.0e-1
sim = test_convergence(dts, prob, KahanLi6(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

sim = test_convergence(dts, prob, McAte8(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 8 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

sim = test_convergence(dts, prob, KahanLi8(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 8 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

dts = 1.0 ./ 2.0 .^ (2:-1:-2)
sim = test_convergence(dts, prob, SofSpa10(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 10 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

################# Out of place symplectic

println("Out of Place")

u0 = 0.0
v0 = 1.0
function f1_harmonic_nip(v, u, p, t)
    return -u
end
function f2_harmonic_nip(v, u, p, t)
    return v
end

ff_harmonic_nip = DynamicalODEFunction(
    f1_harmonic_nip, f2_harmonic_nip;
    analytic = harmonic_analytic
)
prob = DynamicalODEProblem(ff_harmonic_nip, v0, u0, (0.0, 5.0))

sol = solve(prob, SymplecticEuler(), dt = 1 / 10)

dts = 1 .// 2 .^ (6:-1:3)
# Symplectic Euler
sim = test_convergence(dts, prob, SymplecticEuler(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 1 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 1 rtol = 1.0e-1
# Verlet
sim = test_convergence(dts, prob, VelocityVerlet(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1.0e-1
# Test that position converges faster for Verlet
position_error = :final => [
    mean(sim[i].u[2].x[1] - sim[i].u_analytic[2].x[1])
        for i in 1:length(sim)
]
@test first(DiffEqDevTools.calc𝒪estimates(position_error).second) ≈ 4.0 rtol = 1.0e-1

# 2nd Order Tableaus
sim = test_convergence(dts, prob, VerletLeapfrog(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1.0e-1
sim = test_convergence(dts, prob, LeapfrogDriftKickDrift(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1.0e-1
sim = test_convergence(dts, prob, PseudoVerletLeapfrog(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1.0e-1
sim = test_convergence(dts, prob, McAte2(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1.0e-1

# Ruth
sim = test_convergence(dts, prob, Ruth3(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 3 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 3 rtol = 1.0e-1
sim = test_convergence(dts, prob, McAte3(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 3 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 3 rtol = 1.0e-1
sim = test_convergence(dts, prob, CandyRoz4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, McAte4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, McAte42(), dense_errors = true)
@test_broken sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test_broken sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, CalvoSanz4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

dts = 1 .// 2 .^ (4:-1:0)
sim = test_convergence(dts, prob, McAte5(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

sim = test_convergence(dts, prob, Yoshida6(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4.69 rtol = 1.0e-1
sim = test_convergence(dts, prob, KahanLi6(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

sim = test_convergence(dts, prob, McAte8(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 8 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

sim = test_convergence(dts, prob, KahanLi8(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 8 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

dts = 1.0 ./ 2.0 .^ (2:-1:-2)
sim = test_convergence(dts, prob, SofSpa10(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 10 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

################# f1 dependent on v

println("f1 dependent on v")

u0 = fill(0.0, 2)
v0 = ones(2)
function f1_v(dv, v, u, p, t)
    return dv .= v
end
function f2_v(du, v, u, p, t)
    return du .= v
end
function f_v_analytic(y0, p, x)
    v0, u0 = y0.x
    return ArrayPartition(v0 * exp(x), v0 * exp(x) - v0 + u0)
end
ff_v = DynamicalODEFunction(f1_v, f2_v; analytic = f_v_analytic)
prob = DynamicalODEProblem(ff_v, v0, u0, (0.0, 5.0))

dts = 1 .// 2 .^ (6:-1:3)
# LeapfrogDriftKickDrift
sim = test_convergence(dts, prob, LeapfrogDriftKickDrift(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1.0e-1
