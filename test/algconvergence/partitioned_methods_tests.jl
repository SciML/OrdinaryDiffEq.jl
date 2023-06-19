using OrdinaryDiffEq, Test, RecursiveArrayTools, DiffEqDevTools, Statistics

u0 = fill(0.0, 2)
v0 = ones(2)
function f1_harmonic(dv, v, u, p, t)
    dv .= -u
end
function f2_harmonic(du, v, u, p, t)
    du .= v
end
function harmonic_analytic(y0, p, x)
    v0, u0 = y0.x
    ArrayPartition(-u0 * sin(x) + v0 * cos(x), u0 * cos(x) + v0 * sin(x))
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

@test sol[end][1] == sol2[end][1]
@test sol_verlet[end][1] == sol2_verlet[end][1]
@test sol_ruth3[end][1] == sol2_ruth3[end][1]
@test sol[end][3] == sol2[end][3]
@test sol_verlet[end][3] == sol2_verlet[end][3]
@test sol_ruth3[end][3] == sol2_ruth3[end][3]

prob = DynamicalODEProblem(ff_harmonic, v0, u0, (0.0, 5.0))
println("Convergence tests")

dts = 1 .// 2 .^ (6:-1:3)
# Symplectic Euler
sim = test_convergence(dts, prob, SymplecticEuler(), dense_errors = true)
@test sim.𝒪est[:l2]≈1 rtol=1e-1
@test sim.𝒪est[:L2]≈1 rtol=1e-1
# Verlet
sim = test_convergence(dts, prob, VelocityVerlet(), dense_errors = true)
@test sim.𝒪est[:l2]≈2 rtol=1e-1
@test sim.𝒪est[:L2]≈2 rtol=1e-1
# Test that position converges faster for Verlet
position_error = :final => [mean(sim[i].u[2].x[1] - sim[i].u_analytic[2].x[1])
                            for i in 1:length(sim)]
@test first(DiffEqDevTools.calc𝒪estimates(position_error).second)≈4.0 rtol=1e-1

# 2nd Order Tableaus
sim = test_convergence(dts, prob, VerletLeapfrog(), dense_errors = true)
@test sim.𝒪est[:l2]≈2 rtol=1e-1
@test sim.𝒪est[:L2]≈2 rtol=1e-1
sim = test_convergence(dts, prob, PseudoVerletLeapfrog(), dense_errors = true)
@test sim.𝒪est[:l2]≈2 rtol=1e-1
@test sim.𝒪est[:L2]≈2 rtol=1e-1
sim = test_convergence(dts, prob, McAte2(), dense_errors = true)
@test sim.𝒪est[:l2]≈2 rtol=1e-1
@test sim.𝒪est[:L2]≈2 rtol=1e-1

# Ruth
sim = test_convergence(dts, prob, Ruth3(), dense_errors = true)
@test sim.𝒪est[:l2]≈3 rtol=1e-1
@test sim.𝒪est[:L2]≈3 rtol=1e-1
sim = test_convergence(dts, prob, McAte3(), dense_errors = true)
@test sim.𝒪est[:l2]≈3 rtol=1e-1
@test sim.𝒪est[:L2]≈3 rtol=1e-1
sim = test_convergence(dts, prob, CandyRoz4(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, McAte4(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, McAte42(), dense_errors = true)
@test_broken sim.𝒪est[:l2]≈4 rtol=1e-1
@test_broken sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, CalvoSanz4(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

dts = 1 .// 2 .^ (4:-1:0)
sim = test_convergence(dts, prob, McAte5(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

sim = test_convergence(dts, prob, Yoshida6(), dense_errors = true)
@test sim.𝒪est[:l2]≈6 rtol=1e-1
@test sim.𝒪est[:L2]≈4.69 rtol=1e-1
sim = test_convergence(dts, prob, KahanLi6(), dense_errors = true)
@test sim.𝒪est[:l2]≈6 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

sim = test_convergence(dts, prob, McAte8(), dense_errors = true)
@test sim.𝒪est[:l2]≈8 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

sim = test_convergence(dts, prob, KahanLi8(), dense_errors = true)
@test sim.𝒪est[:l2]≈8 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

sol = solve(prob, Nystrom4(), dt = 1 / 1000)

# Nyström method
dts = 1 .// 2 .^ (9:-1:6)
sim = test_convergence(dts, prob, Nystrom4(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, Nystrom4VelocityIndependent(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, IRKN3(), dense_errors = true)
@test sim.𝒪est[:l2]≈3 rtol=1e-1
@test sim.𝒪est[:L2]≈3 rtol=1e-1
sim = test_convergence(dts, prob, IRKN4(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
dts = 1.0 ./ 2.0 .^ (5:-1:0)
sim = test_convergence(dts, prob, Nystrom5VelocityIndependent(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, FineRKN4(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, FineRKN5(), dense_errors = true)
@test sim.𝒪est[:l2]≈6 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

dts = 1.0 ./ 2.0 .^ (2:-1:-2)
sim = test_convergence(dts, prob, SofSpa10(), dense_errors = true)
@test sim.𝒪est[:l2]≈10 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

# Methods need BigFloat to test convergence rate
dts = big"1.0" ./ big"2.0" .^ (5:-1:1)
prob_big = DynamicalODEProblem(ff_harmonic, [big"1.0", big"1.0"],
    [big"0.0", big"0.0"], (big"0.", big"70."))
sim = test_convergence(dts, prob_big, DPRKN4(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, DPRKN5(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, DPRKN6(), dense_errors = true)
@test sim.𝒪est[:l2]≈6 rtol=1e-1
@test sim.𝒪est[:L2]≈6 rtol=1e-1
sim = test_convergence(dts, prob_big, DPRKN6FM(), dense_errors = true)
@test sim.𝒪est[:l2]≈6 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, DPRKN8(), dense_errors = true)
@test sim.𝒪est[:l2]≈8 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, DPRKN12(), dense_errors = true)
@test sim.𝒪est[:l2]≈12 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, ERKN4(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, ERKN5(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, ERKN7(), dense_errors = true)
@test sim.𝒪est[:l2]≈7 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

# Adaptive methods regression test
sol = solve(prob, FineRKN4())
@test length(sol.u) < 16
sol = solve(prob, DPRKN4())
@test length(sol.u) < 25
sol = solve(prob, DPRKN5())
@test length(sol.u) < 38
sol = solve(prob, DPRKN6())
@test length(sol.u) < 20
sol = solve(prob, DPRKN6FM())
@test length(sol.u) < 25
sol = solve(prob, DPRKN8())
@test length(sol.u) < 13
sol = solve(prob, DPRKN12())
@test length(sol.u) < 10
sol = solve(prob, ERKN4(), reltol = 1e-8)
@test length(sol.u) < 38
sol = solve(prob, ERKN5(), reltol = 1e-8)
@test length(sol.u) < 34
sol = solve(prob, ERKN7(), reltol = 1e-8)
@test length(sol.u) < 38

# Test array partition outside of symplectic

function f_ap(du, u, p, t)
    du.x[1] .= -2u.x[2]
    du.x[2] .= u.x[1]
end

u = ArrayPartition((u0, v0))

prob = ODEProblem(f_ap, u, (0.0, 5.0))
sol = solve(prob, Euler(), dt = 1 / 100)

################# Out of place symplectic

println("Out of Place")

u0 = 0.0
v0 = 1.0
function f1_harmonic_nip(v, u, p, t)
    -u
end
function f2_harmonic_nip(v, u, p, t)
    v
end

ff_harmonic_nip = DynamicalODEFunction(f1_harmonic_nip, f2_harmonic_nip;
    analytic = harmonic_analytic)
prob = DynamicalODEProblem(ff_harmonic_nip, v0, u0, (0.0, 5.0))

sol = solve(prob, SymplecticEuler(), dt = 1 / 10)

dts = 1 .// 2 .^ (6:-1:3)
# Symplectic Euler
sim = test_convergence(dts, prob, SymplecticEuler(), dense_errors = true)
@test sim.𝒪est[:l2]≈1 rtol=1e-1
@test sim.𝒪est[:L2]≈1 rtol=1e-1
# Verlet
sim = test_convergence(dts, prob, VelocityVerlet(), dense_errors = true)
@test sim.𝒪est[:l2]≈2 rtol=1e-1
@test sim.𝒪est[:L2]≈2 rtol=1e-1
# Test that position converges faster for Verlet
position_error = :final => [mean(sim[i].u[2].x[1] - sim[i].u_analytic[2].x[1])
                            for i in 1:length(sim)]
@test first(DiffEqDevTools.calc𝒪estimates(position_error).second)≈4.0 rtol=1e-1

# 2nd Order Tableaus
sim = test_convergence(dts, prob, VerletLeapfrog(), dense_errors = true)
@test sim.𝒪est[:l2]≈2 rtol=1e-1
@test sim.𝒪est[:L2]≈2 rtol=1e-1
sim = test_convergence(dts, prob, PseudoVerletLeapfrog(), dense_errors = true)
@test sim.𝒪est[:l2]≈2 rtol=1e-1
@test sim.𝒪est[:L2]≈2 rtol=1e-1
sim = test_convergence(dts, prob, McAte2(), dense_errors = true)
@test sim.𝒪est[:l2]≈2 rtol=1e-1
@test sim.𝒪est[:L2]≈2 rtol=1e-1

# Ruth
sim = test_convergence(dts, prob, Ruth3(), dense_errors = true)
@test sim.𝒪est[:l2]≈3 rtol=1e-1
@test sim.𝒪est[:L2]≈3 rtol=1e-1
sim = test_convergence(dts, prob, McAte3(), dense_errors = true)
@test sim.𝒪est[:l2]≈3 rtol=1e-1
@test sim.𝒪est[:L2]≈3 rtol=1e-1
sim = test_convergence(dts, prob, CandyRoz4(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, McAte4(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, McAte42(), dense_errors = true)
@test_broken sim.𝒪est[:l2]≈4 rtol=1e-1
@test_broken sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, CalvoSanz4(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

dts = 1 .// 2 .^ (4:-1:0)
sim = test_convergence(dts, prob, McAte5(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

sim = test_convergence(dts, prob, Yoshida6(), dense_errors = true)
@test sim.𝒪est[:l2]≈6 rtol=1e-1
@test sim.𝒪est[:L2]≈4.69 rtol=1e-1
sim = test_convergence(dts, prob, KahanLi6(), dense_errors = true)
@test sim.𝒪est[:l2]≈6 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

sim = test_convergence(dts, prob, McAte8(), dense_errors = true)
@test sim.𝒪est[:l2]≈8 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

sim = test_convergence(dts, prob, KahanLi8(), dense_errors = true)
@test sim.𝒪est[:l2]≈8 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

# Nyström method
dts = 1 .// 2 .^ (9:-1:6)
sim = test_convergence(dts, prob, Nystrom4(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, Nystrom4VelocityIndependent(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
@test_broken sim = test_convergence(dts, prob, IRKN3(), dense_errors = true)
@test_broken sim.𝒪est[:l2]≈3 rtol=1e-1
@test_broken sim.𝒪est[:L2]≈3 rtol=1e-1
@test_broken sim = test_convergence(dts, prob, IRKN4(), dense_errors = true)
#@test_broken sim.𝒪est[:l2] ≈ 4 rtol = 1e-1
#@test_broken sim.𝒪est[:L2] ≈ 4 rtol = 1e-1
dts = 1.0 ./ 2.0 .^ (5:-1:0)
sim = test_convergence(dts, prob, Nystrom5VelocityIndependent(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, FineRKN4(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob, FineRKN5(), dense_errors = true)
@test sim.𝒪est[:l2]≈6 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

dts = 1.0 ./ 2.0 .^ (2:-1:-2)
sim = test_convergence(dts, prob, SofSpa10(), dense_errors = true)
@test sim.𝒪est[:l2]≈10 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

# Methods need BigFloat to test convergence rate
dts = big"1.0" ./ big"2.0" .^ (5:-1:1)
prob_big = DynamicalODEProblem(ff_harmonic_nip, big"1.0", big"0.0",
    (big"0.", big"70."))
sim = test_convergence(dts, prob_big, DPRKN4(), dense_errors = true)
@test sim.𝒪est[:l2]≈4 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, DPRKN5(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, DPRKN6(), dense_errors = true)
@test sim.𝒪est[:l2]≈6 rtol=1e-1
@test sim.𝒪est[:L2]≈6 rtol=1e-1
sim = test_convergence(dts, prob_big, DPRKN6FM(), dense_errors = true)
@test sim.𝒪est[:l2]≈6 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, DPRKN8(), dense_errors = true)
@test sim.𝒪est[:l2]≈8 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, DPRKN12(), dense_errors = true)
@test sim.𝒪est[:l2]≈12 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, ERKN4(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, ERKN5(), dense_errors = true)
@test sim.𝒪est[:l2]≈5 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1
sim = test_convergence(dts, prob_big, ERKN7(), dense_errors = true)
@test sim.𝒪est[:l2]≈7 rtol=1e-1
@test sim.𝒪est[:L2]≈4 rtol=1e-1

# Adaptive methods regression test
sol = solve(prob, FineRKN4())
@test length(sol.u) < 16
sol = solve(prob, DPRKN4())
@test length(sol.u) < 25
sol = solve(prob, DPRKN5())
@test length(sol.u) < 38
sol = solve(prob, DPRKN6())
@test length(sol.u) < 20
sol = solve(prob, DPRKN6FM())
@test length(sol.u) < 25
sol = solve(prob, DPRKN8())
@test length(sol.u) < 13
sol = solve(prob, DPRKN12())
@test length(sol.u) < 10
sol = solve(prob, ERKN4(), reltol = 1e-8)
@test length(sol.u) < 38
sol = solve(prob, ERKN5(), reltol = 1e-8)
@test length(sol.u) < 34
sol = solve(prob, ERKN7(), reltol = 1e-8)
@test length(sol.u) < 38
