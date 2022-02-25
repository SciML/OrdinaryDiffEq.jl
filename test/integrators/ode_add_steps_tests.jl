using OrdinaryDiffEq, Test

function test_ode(u, p, t)
    [p[1] - (1 - p[1])*u[1]]
end

function test_ode(du, u, p, t)
    du[1] = p[1] - (1 - p[1])*u[1]
    return nothing
end

test_solution(t) = t <= 5 ? t : 5. * exp(-(t-5))

tspan = (0.,10.)
testtimes = range(tspan[1], stop=tspan[2], length=1001)
pullback_condition(u, t, i) = t - 5
pullback_affect!(i) = i.p[1] = abs(1 - i.p[1])
cb = ContinuousCallback(pullback_condition, pullback_affect!)

# DPRKN6 is left out because second order ODE
algs = [Tsit5, Rosenbrock23, DP5, DP8, SSPRK43, SSPRK432, OwrenZen3, SSPRK22, SSPRK33,
        OwrenZen4, OwrenZen5, Rosenbrock32, Rodas5, Rodas4, Rodas42, Rodas5P] ## Works for these
lazy_alg = [BS5, Vern6, Vern7, Vern8, Vern9]

nonstandard_interp_algs = union(algs,bad_algs,lazy_alg)

passed = fill(false, 2length(nonstandard_interp_algs))

cur_itr = 0
for inplace in [false,true], alg in algs
    global cur_itr
    prob = ODEProblem{inplace}(test_ode, [0.], tspan, [1.])
    sol = solve(prob, alg(); callback=cb,dt=0.0013)
    pass = all(isapprox(sol(t)[1], test_solution(t); atol=0.05) for t in testtimes)
    x = maximum(sol(t)[1] - test_solution(t) for t in testtimes)
    cur_itr += 1
    passed[cur_itr] = pass
end

for inplace in [false,true], alg in lazy_alg
    global cur_itr
    prob = ODEProblem{inplace}(test_ode, [0.], tspan, [1.])
    sol = solve(prob, alg(); callback=cb,dt=0.0013)
    fail = all(isapprox(sol(t)[1], test_solution(t); atol=0.05) for t in testtimes)

    prob = ODEProblem{inplace}(test_ode, [0.], tspan, [1.])
    sol = solve(prob, alg(lazy=false); callback=cb,dt=0.0013)
    pass = all(isapprox(sol(t)[1], test_solution(t); atol=0.05) for t in testtimes)

    cur_itr += 1
    passed[cur_itr] = pass && !fail
end

any(.!(passed)) && @warn("The following algorithms failed the continuous callback test: $(vcat(algs,algs,lazy_alg,lazy_alg)[.!(passed)])")
@test all(passed)
